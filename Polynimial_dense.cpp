#include <algorithm>
#include <string>
#include <iostream>
#include<vector>
#include<set>
#include <iterator>
using namespace std;

template <typename T>
class Polynomial {
private:
	vector<T> coef;
public:
	Polynomial(const vector<T>& coefficients) {
		coef = coefficients;
		Cutting();
	}
	Polynomial(const T value = T()) {
		coef.push_back(value);
		Cutting();
	}
	template <class It>
	Polynomial(It first, It last) {
		T next_coef;
		while (first != last) {
			next_coef = *first;
			coef.push_back(next_coef);
			++first;
		}
		Cutting();
	}

	vector<T> get_coef() const {
		return coef;
	}
	void print_coef() const {
		for (auto el : coef) {
			cout << el << "  ";
		}
		cout << endl;
	}
	void Cutting() {
		while (!coef.empty() && coef.back() == T(0)) {
			coef.pop_back();
		}
	}
	int Degree() const {
		return (coef.size() - 1);
	}
	void normalization() {
		for (int i = 0; i <= this->Degree(); ++i) {
			coef[i] /= coef[this->Degree()];
		}
	}
	T module(T num) {
		if (num > T(0)) {
			return num;
		}
		else {
			return num * T(-1);
		}
	}
	Polynomial<T>& operator += (const Polynomial<T>& p) {
		int new_degree = max(this->Degree(), p.Degree());
		coef.resize(new_degree + 1);
		for (int i = 0; i <= p.Degree(); ++i) {
			auto tmp = p[i];
			coef[i] += tmp;
		}
		Cutting();
		return *this;
	}
	Polynomial<T>& operator -= (const Polynomial<T>& p) {
		int new_degree = max(this->Degree(), p.Degree());
		coef.resize(new_degree + 1);
		for (int i = 0; i < p.Degree() + 1; ++i) {
			coef[i] -= p.get_coef()[i];
		}
		Cutting();
		return *this;
	}
	Polynomial<T>& operator *= (const Polynomial& p) {
		vector<T> result_coef = vector<T>(this->Degree() + p.Degree() + 1);
		for (int i = 0; i <= this->Degree(); ++i) {
			for (int j = 0; j <= p.Degree(); ++j) {
				result_coef[i + j] += coef[i] * p[j];
			}
		}
		Polynomial<T> result = Polynomial<T>(result_coef);
		result.Cutting();
		*this = result;
		return *this;
	}

	typename vector<T>::const_iterator begin() const {
		return coef.begin();
	}
	typename vector<T>::const_iterator end() const {
		return coef.end();
	}

	T operator [] (int pos) const {
		if (pos > this->Degree()) {
			return T(0);
		}
		else {
			return coef[pos];
		}
	}
	T operator() (T x) const {
		T result = T(0);
		int ind = this->Degree();
		while (ind >= 0) {
			result = result*x + coef[ind];
			--ind;
		}
		return result;
	}

	bool operator == (const Polynomial<T>& p) const {
		return coef == p.get_coef();
	}
	bool operator != (const Polynomial<T>& p) const {
		return coef != p.get_coef();
	}

	friend Polynomial<T> operator + (Polynomial<T> p1, const Polynomial<T>& p2) {
		return p1 += p2;
	}
	friend Polynomial<T> operator - (Polynomial<T> p1, const Polynomial<T>& p2) {
		return p1 -= p2;
	}
	friend Polynomial<T> operator * (Polynomial<T> p1, const Polynomial<T>& p2) {
		return p1 *= p2;
	}

	friend Polynomial<T> operator & (const Polynomial<T> p1, const Polynomial<T>& p2) {
		Polynomial<T> tmp = Polynomial<T>(T(0));
		Polynomial<T> composition;
		for (int i = 0; i <= p1.Degree(); ++i) {
			tmp = Polynomial<T>(p1[i]);
			for (int j = 0; j < i; ++j) {
				tmp *= p2;
			}
			composition += tmp;
		}
		return composition;
	}

	friend Polynomial<T> operator / (const Polynomial<T>& p1, const Polynomial<T>& p2) {
		if (p2.Degree() != -1) {
			Polynomial<T> p1_copy = p1;	
			int quotient_size = max(0, p1_copy.Degree() - p2.Degree() + 1);
			vector<T> quotient_coef(quotient_size);
			vector<T> monom(quotient_size, T(0));
			Polynomial<T> monomial = Polynomial<T>(monom);
			T monom_coef;
			int cur_quotient_ind;
			while (p1_copy.Degree() >= p2.Degree())
			{
				monom_coef = p1_copy[p1_copy.Degree()] / p2[p2.Degree()];
				cur_quotient_ind = p1_copy.Degree() - p2.Degree();
				quotient_coef[cur_quotient_ind] = monom_coef;
				monom[cur_quotient_ind] = monom_coef;
				monomial = Polynomial<T>(monom);
				p1_copy -= (monomial*p2);
				monom[cur_quotient_ind] = T(0);
			}
			return Polynomial<T>(quotient_coef);
		}
		else {
			return{};
		}
	}
	friend Polynomial<T> operator % (const Polynomial<T>& p1, const Polynomial<T>& p2) {
		return p1 - (p1 / p2)*p2;
	}
	friend Polynomial<T> operator , (const Polynomial<T>& p1, const Polynomial<T>& p2) {
		Polynomial<T> p1_copy = p1;
		Polynomial<T> p2_copy = p2;
		while (p1_copy.Degree() != -1 && p2_copy.Degree() != -1) {
			if (p1_copy.Degree() > p2_copy.Degree()) {
				p1_copy = p1_copy % p2_copy;
			}
			else if (p1_copy.Degree() < p2_copy.Degree()) {
				p2_copy = p2_copy % p1_copy;
			}
			else {
				if (p1_copy.module(p1_copy[p1_copy.Degree()]) > p2_copy.module(p2_copy[p2_copy.Degree()])) {
					p1_copy = p1_copy % p2_copy;
				}
				else {
					p2_copy = p2_copy % p1_copy;
				}
			}
		}
		if (p1_copy.Degree() == -1) {
			p2_copy.normalization();
			return p2_copy;
		}
		p1_copy.normalization();
		return p1_copy;
	}
};

template <typename T>
std::ostream& operator << (std::ostream& out, const Polynomial<T>& p) {
	Polynomial<T> p_copy = p;
	p_copy.Cutting();
	if (p_copy.Degree() == -1) {
		out << T(0);
	}
	else {
		string plus_str;
		for (int i = p.Degree(); i >= 0; --i) {
			if (p[i] != T(0)) {
				if (p[i] > T(0)) {
					out << plus_str;;
				}
				if (i > 1) {
					if (p[i] * p[i] != T(1)) {
						out << p[i] << "*";
					}
					else if (p[i] == T(-1)) {
						out << "-";
					}
					out << "x^" << i;
				}
				else {
					if (i == 1) {
						if (p[i] * p[i] != T(1)) {
							out << p[i] << "*";
						}
						else if (p[i] == T(-1)) {
							out << "-";
						}
						out << "x";
					}
					else {
						out << p[i];
					}
				}
			}
			plus_str = "+";
		}
	}
	return out;
}

int main()
{
    return 0;
}
