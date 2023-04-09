#ifndef BIGINTEGER_BIGINTEGER_H
#define BIGINTEGER_BIGINTEGER_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

static const int MOD = 1000000000;

class BigInteger {
public:
    BigInteger(): digits({0}), is_positive(true) {}
    BigInteger (int number) {
        is_positive = number >= 0;
        if (!is_positive) {
            number *= -1;
        }
        digits.push_back(number % MOD);
        if (number >= MOD)
            digits.push_back(number / MOD);
        delete_zeros();
    }
    explicit BigInteger(const std::string& number) {
        int is_negative = number[0] == '-';
        is_positive = !is_negative;
        for (int i = static_cast<int> (number.size()); i > is_negative; i -= 9) {
            if (i < 9) {
                digits.push_back(atoi(number.substr(is_negative, i - is_negative).c_str()));
            } else {
                digits.push_back(atoi(number.substr(i - 9, 9).c_str()));
            }
        }
    }
    explicit BigInteger(std::vector<int>& digits, bool is_positive = true)
        : digits(digits), is_positive(is_positive) {};
    BigInteger (const BigInteger& other) = default;
    void swap(BigInteger& first, BigInteger& second) {
        std::swap(first.digits, second.digits);
        std::swap(first.is_positive, second.is_positive);
    }
    BigInteger& operator= (const BigInteger& other) = default;
    std::vector<int> get_const_digits() {
        return digits;
    }
    size_t get_size() const {
        return digits.size();
    }
    bool get_is_positive() const {
        return is_positive;
    }
    std::string toString () const {
        std::string result;
        for (size_t i = 0; i < digits.size(); ++i) {
            int ten_degree = 1;
            for (int j = 0; j < 9; ++j) {
                result.push_back(static_cast<char>(((digits[i] / ten_degree) % 10) + '0') );
                ten_degree *= 10;
            }
        }
        while (result.size() > 1 && result.back() == '0') {
            result.pop_back();
        }
        std::string answer;
        if (!is_positive) {
            answer.push_back('-');
        }
        for (size_t i = 0; i < result.size(); ++i) {
            answer.push_back(result[result.size() - i - 1]);
        }
        return answer;
    }
    void set_is_positive (bool value) {
        is_positive = value;
    }
    bool operator < (const BigInteger& other) const {
        if (is_positive != other.get_is_positive()) {
            return !is_positive;
        }
        if (digits.size() < other.get_size()) {
            return is_positive;
        } else if (digits.size() > other.get_size()) {
            return !is_positive;
        }
        for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
            if (digits[i] > other.digits[i]) {
                return !is_positive;
            } else if (digits[i] < other.digits[i]) {
                return is_positive;
            }
        }
        return false;
    }
    bool operator > (const BigInteger& other) const {
        return (other < (*this));
    }
    bool operator >= (const BigInteger& other) const {
        return !((*this) < other);
    }
    bool operator <= (const BigInteger& other) const {
        return !(other < (*this));
    }
    bool operator == (const BigInteger& other) const {
        return (!((*this) < other) && !(other < (*this)));
    }
    bool operator != (const BigInteger& other) const {
        return ((*this) < other || other < (*this));
    }
    explicit operator bool() {
        return !(digits.size() == 1 && digits[0] == 0);
    }

    BigInteger& operator++ () {
        if (is_positive)
            abs_sum(BigInteger(0), true);
        else
            abs_diff(BigInteger(0), true);
        return *this;
    }
    BigInteger operator++ (int) {
        BigInteger temporary = *this;
        ++(*this);
        return temporary;
    }
    BigInteger& operator-- () {
        if (is_positive)
            abs_diff(BigInteger(0), true);
        else
            abs_sum(BigInteger(0), true);
        return *this;
    }
    BigInteger operator-- (int) {
        BigInteger temporary = *this;
        --(*this);
        return temporary;
    }
    BigInteger operator - () {
        BigInteger tmp = *this;
        if (tmp.digits.size() > 1 || tmp.digits[0] != 0) {
            tmp.is_positive = !tmp.is_positive;
        }
        return tmp;
    }

    BigInteger& operator*= (int number) {
        size_t size = digits.size();
        for (size_t i = 0, carry = 0; i < size || carry; ++i) {
            if (digits.size() == i) {
                digits.push_back(0);
            }
            long long res = digits[i] * 1LL * number + carry;
            digits[i] = static_cast<int> (res % MOD);
            carry = static_cast<int> (res / MOD);
        }
        return (*this);
    }
    BigInteger& operator+= (const BigInteger& other) {
        if (is_positive && !other.is_positive) {
            is_positive = abs_diff(other, false);
        } else if (!is_positive && other.get_is_positive()) {
            is_positive = !abs_diff(other, false);
        } else {
            abs_sum(other, false);
        }
        return (*this);
    }
    BigInteger& operator-= (const BigInteger& other) {
        if (is_positive != other.get_is_positive()) {
            abs_sum(other, false);
        } else if (is_positive) {
            is_positive = abs_diff(other, false);
        } else {
            is_positive = !abs_diff(other, false);
        }
        return (*this);
    }

    BigInteger& operator *= (const BigInteger& other) {
        if (!other.is_positive) {
            is_positive = !is_positive;
        }
        size_t old_size = digits.size();
        size_t old_other_size = other.digits.size();
        for (size_t i = 0; i < old_other_size; ++i) {
            digits.push_back(0);
        }
        for (int i = static_cast<int>(old_size) - 1; i >= 0; --i) {
            int addition = 0;
            int digit = digits[i];
            for (size_t j = 0; j < other.digits.size() || addition; ++j) {
                long long res = digit * 1LL * (j < other.get_size() ? other.digits[j] : 0)
                                + addition + (j != 0 ? digits[i+j] : 0);
                digits[i + j] = static_cast<int> (res % MOD);
                addition = static_cast<int> (res / MOD);
            }
        }
        delete_zeros();
        if (digits.size() == 1 && digits[0] == 0) {
            is_positive = true;
        }
        return (*this);
    }
    BigInteger& operator%=(const BigInteger &other) {
        for (int i = static_cast<int>(digits.size()) - other.get_size(); i >= 0; --i) {
            if (!abs_less_shifted(other, i)) {
                int lhs = 0;
                int rhs = MOD;
                while (rhs - lhs > 1) {
                    int mid = (lhs + rhs) >> 1;
                    if (!abs_less_shifted(other * mid, i)) {
                        lhs = mid;
                    } else {
                        rhs = mid;
                    }
                }
                difference (other * lhs, i);
            }
        }
        delete_zeros();
        if (digits.size() == 1 && digits[0] == 0) {
            is_positive = true;
        }
        return *this;
    }
    BigInteger& operator/= (const BigInteger &other) {
        BigInteger result;
        is_positive = !(is_positive ^ other.is_positive);
        for (int i = static_cast<int>(digits.size()) - other.get_size(); i >= 0; --i) {
            int lhs = 0;
            if (!abs_less_shifted(other, i)) {
                int rhs = MOD;
                while (rhs - lhs > 1) {
                    int mid = (lhs + rhs) >> 1;
                    if (!abs_less_shifted(other * mid, i)) {
                        lhs = mid;
                    } else {
                        rhs = mid;
                    }
                }
                difference(other * lhs, i);
            }
            result.digits.push_back(lhs);
        }
        std::reverse(result.digits.begin(), result.digits.end());
        std::swap(digits, result.digits);
        delete_zeros();
        if (digits.size() == 1 && digits[0] == 0) {
            is_positive = true;
        }
        return *this;
    }

    BigInteger operator* (int other) const {
        BigInteger tmp(*this);
        tmp *= BigInteger(other);
        return tmp;
    }
    void delete_zeros () {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
    }
    bool is_even() {
        return (digits[0] % 2) == 0;
    }
    bool is_equal_to (int num) {
        return (digits.size() == 1 && digits[0] == num);
    }
private:
    std::vector<int> digits;
    bool is_positive;
    bool abs_compare_less (const BigInteger& other) const {
        if (digits.size() < other.get_size()) {
            return true;
        }
        if (digits.size() > other.get_size()) {
            return false;
        }
        for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
            if (digits[i] > other.digits[i]) {
                return false;
            } else if (digits[i] < other.digits[i]) {
                return true;
            }
        }
        return false;
    }
    void difference (const BigInteger& other, int digits_to_shift) {
        int addition = 0;
        int digit;
        size_t size = digits.size();
        for (size_t i = digits_to_shift; i < size || addition; ++i) {
            if (i == digits.size()) {
                digits.push_back(0);
            }
            digit = digits[i] - addition - (i - digits_to_shift < other.get_size()
                                            ? other.digits[i - digits_to_shift] : 0);
            if (digit < 0) {
                digits[i] = MOD + digit;
                addition = 1;
            } else {
                digits[i] = digit;
                addition = 0;
            }
        }
        delete_zeros();
    }
    bool abs_less_shifted (const BigInteger& other, int digits_to_shift = 0) {
        if (get_size() < other.get_size() + digits_to_shift) {
            return true;
        } else if (get_size() > other.get_size() + digits_to_shift) {
            return false;
        }
        for (int i = static_cast<int> (other.get_size()) + digits_to_shift - 1; i >= digits_to_shift; --i) {
            if (digits[i] > other.digits[i-digits_to_shift]) {
                return false;
            } else if (digits[i] < other.digits[i-digits_to_shift]) {
                return true;
            }
        }
        return false;
    }
    void abs_sum (const BigInteger& other, bool is_increment) {
        int addition = is_increment;
        size_t max_size = std::max(digits.size(), other.digits.size());
        for (size_t i = 0; (i < max_size && !is_increment) || addition; ++i) {
            if (i == digits.size())
                digits.push_back (0);
            digits[i] += addition + (i < other.digits.size() ? other.digits[i] : 0);
            addition = digits[i] >= MOD;
            if (addition)  digits[i] -= MOD;
        }
    }
    // return whether result is positive
    bool abs_diff (const BigInteger& other, bool is_increment) {
        bool is_answer_negative = abs_compare_less(other);
        int addition = is_increment;
        int digit;
        size_t max_size = std::max(digits.size(), other.digits.size());
        for (size_t i = 0; (i < max_size && !is_increment) || addition; ++i) {
            if (digits.size() == i) {
                digits.push_back(0);
            }
            if (is_answer_negative) {
                digit = other.digits[i] - digits[i] - addition;
            } else {
                digit = digits[i] - addition - (i < other.get_size() ? other.digits[i] : 0);
            }
            if (digit < 0) {
                digits[i] = MOD + digit;
                addition = 1;
            } else {
                digits[i] = digit;
                addition = 0;
            }
        }
        delete_zeros();
        return !is_answer_negative;
    }
};

BigInteger operator+ (const BigInteger& left, const BigInteger& right) {
    BigInteger result (left);
    return (result += right);
}
BigInteger operator- (const BigInteger& left, const BigInteger& right) {
    BigInteger result (left);
    return (result -= right);
}
BigInteger operator* (const BigInteger& left, const BigInteger& right) {
    BigInteger result (left);
    return (result *= right);
}
BigInteger operator% (const BigInteger& left, const BigInteger& right) {
    BigInteger result (left);
    return (result %= right);
}
BigInteger operator/ (const BigInteger& left, const BigInteger& right) {
    BigInteger result (left);
    return (result /= right);
}
std::istream& operator >> (std::istream& in, BigInteger& number) {
    std::string input;
    in >> input;
    number = BigInteger(input);
    return in;
}
std::ostream& operator << (std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

class Rational {
public:
    Rational() = default;
    Rational(const BigInteger& numerator_, const BigInteger& denominator_)
            : numerator(numerator_),  denominator(denominator_) {
        make_positive();
        coprime();
    }
    Rational (const BigInteger& numerator):numerator(numerator) {
        denominator = BigInteger(1);
    }
    explicit Rational (int numerator, int denominator)
            : numerator(BigInteger(numerator)), denominator(BigInteger(denominator)) {
        make_positive();
        coprime();
    }
    Rational(int numerator): numerator(numerator), denominator(1) {};
    Rational& operator+= (const Rational& other) {
        numerator = numerator * other.denominator
                                   + denominator * other.numerator;
        denominator = denominator * other.denominator;
        make_positive();
        coprime();
        return *this;
    }
    Rational& operator-= (const Rational& other) {
        numerator = numerator * other.denominator
                                   - denominator * other.numerator;
        denominator = denominator * other.denominator;
        make_positive();
        coprime();
        return (*this);
    }
    Rational& operator*= (const Rational& other) {
        numerator *= other.numerator;
        denominator *= other.denominator;
        make_positive();
        coprime();
        return (*this);
    }
    Rational& operator/= (const Rational& other) {
        numerator *= other.denominator;
        denominator *= other.numerator;
        make_positive();
        coprime();
        return *this;
    }
    bool operator < (const Rational& other) const {
        return (numerator * other.denominator < other.numerator * denominator);
    }
    bool operator > (const Rational& other) const {
        return (other < (*this));
    }
    bool operator <= (const Rational& other) const {
        return !(*this < other);
    }
    bool operator >= (const Rational& other) const {
        return !(other < *this);
    }
    // may be much faster
    bool operator == (const Rational& other) const {
        return other.numerator == numerator && other.denominator == denominator;
    }
    bool operator != (const Rational& other) const {
        return (*this < other || other < *this);
    }
    Rational& operator= (const Rational& other) {
        numerator = other.numerator;
        denominator = other.denominator;
        return (*this);
    }
    // unary minus
    Rational operator- () {
        Rational tmp (numerator);
        if (numerator != 0)
            tmp.numerator.set_is_positive(!numerator.get_is_positive());
        return tmp;
    }
    std::string toString() {
        std::string result;
        if (denominator == 1 || denominator == -1) {
            if (denominator.get_is_positive()) {
                result = numerator.toString();
            } else {
                if (numerator.get_is_positive()) {
                    result = '-' + numerator.toString();
                } else {
                    result = numerator.toString();
                    result = result.substr(1, result.size() - 1);
                }
            }
        } else {
            if (numerator.get_is_positive() ^ denominator.get_is_positive()) {
                result = "-";
            }
            result += numerator.toString() + '/';
            std::string tmp = denominator.toString();
            if (tmp[0] == '-') {
                tmp = tmp.substr(1, tmp.size() - 1);
            }
            result += tmp;
        }
        return result;
    }
    std::string asDecimal (size_t precision_digits = 0) const {
        size_t precision = precision_digits / 9;
        long long ten_degree = 1;
        for (size_t i = 0; i < precision_digits % 9; ++i) {
            ten_degree *= 10;
        }
        std::string answer;
        BigInteger quotient = numerator / denominator;
        quotient.delete_zeros();
        BigInteger remainder = numerator - denominator * quotient;
        remainder.delete_zeros();
        if (numerator.get_is_positive() ^ denominator.get_is_positive()) {
            answer = '-';
        }
        std::vector<int> new_numerator_digits = remainder.get_const_digits();
        new_numerator_digits.resize(remainder.get_size() + precision);
        for (int i = numerator.get_size() + precision - 1; i >= 0; --i) {
            if (i >= static_cast<int>(precision)) {
                new_numerator_digits[i] = new_numerator_digits[i-static_cast<int>(precision)];
            } else {
                new_numerator_digits[i] = 0;
            }
        }
        long long tail = ten_degree * new_numerator_digits.back();
        new_numerator_digits.back() = tail % MOD;
        if (tail / MOD != 0) {
            new_numerator_digits.push_back(tail / MOD);
        }
        BigInteger new_numerator(new_numerator_digits);
        new_numerator.delete_zeros();
        BigInteger result = new_numerator / denominator;
        result.delete_zeros();
        result.set_is_positive(true);
        std::string after_point = result.toString();
        std::string new_after_point;
        if (precision_digits > after_point.size()) {
            for (size_t i = 0; i < precision_digits - after_point.size(); ++i) {
                 new_after_point += '0';
            }
        }
        new_after_point += after_point;
        new_after_point = new_after_point.substr(0, precision_digits);
        answer += quotient.toString() + '.' + new_after_point;
        return answer;
    }
    explicit operator double() const{
        size_t precision_needed = (denominator.get_size() - numerator.get_size() + 1) * 9;
        return std::stod(asDecimal(std::max(precision_needed, static_cast<size_t>(16)) ));
    }
private:
    BigInteger numerator;
    BigInteger denominator;
    void make_positive() {
        if ((numerator.get_is_positive() || denominator.get_is_positive()) == 0) {
            numerator.set_is_positive(true);
            denominator.set_is_positive(true);
        }
    }
    void coprime () {
        BigInteger tmp_numerator = numerator;
        tmp_numerator.set_is_positive(true);
        BigInteger tmp_denominator = denominator;
        tmp_denominator.set_is_positive(true);
        BigInteger answer(1);
        BigInteger gcd = greatest_common_divisor(tmp_numerator, tmp_denominator, answer);
        gcd.delete_zeros();
        numerator /= gcd;
        denominator /= gcd;
    }

    BigInteger greatest_common_divisor (BigInteger& left, BigInteger& right,
                                        BigInteger& answer) {
        if (left.is_equal_to(0)) {
            return right * answer;
        }
        if (right.is_equal_to(0)) {
            return left * answer;
        }
        if (left == right) {
            return left * answer;
        }
        while (left.is_even() && right.is_even()) {
            left /= BigInteger(2);
            right /= BigInteger(2);
            answer *= BigInteger(2);
        }
        while (left.is_even()) {
            left /= BigInteger(2);
        }

        while (right.is_even()) {
            right /= BigInteger(2);
        }
        if (left < right) {
            left.swap(left, right);
        }
        left -= right;
        left /= BigInteger(2);
        return greatest_common_divisor(right, left, answer);
    }


};

Rational operator* (const Rational& left, const Rational& right) {
    Rational temporary(left);
    temporary *= right;
    return temporary;
}
Rational operator+ (const Rational& left, const Rational& right) {
    Rational temporary(left);
    temporary += right;
    return temporary;
}
Rational operator- (const Rational& left, const Rational& right) {
    Rational temporary(left);
    temporary -= right;
    return temporary;
}
Rational operator/ (const Rational& left, const Rational& right) {
    Rational temporary(left);
    temporary /= right;
    return temporary;
}

#endif //BIGINTEGER_BIGINTEGER_H
