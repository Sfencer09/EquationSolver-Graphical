/*
 Note: unless otherwise stated, the classes in this project have natural orderings that are inconsistent with equals.
 */
package equation.solver.graphical;

final class my_math_utils {

    static final private int rational_limit = 10000; //caution: increasing this value will drastically increase the time to test if a number is rational

    private my_math_utils() {

    }

    static public boolean is_square(double num) {
        return (((int) Math.pow(((int) Math.sqrt(num)), 2)) == num);
    }

    static public boolean is_cube(double num) {
        return (((int) Math.pow(((int) Math.cbrt(num)), 3)) == num);
    }

    static public boolean is_even(double num) {
        return (num % 2 == 0);
    }

    static public boolean is_negative(double num) {
        return (num < 0);
    }

    static public boolean is_integer(double num) {
        return ((long) num) == num;
    }

    static public int rational_base(double num, int limit) {
        if (Double.isInfinite(num) || Double.isNaN(num)) {
            return 0;
        }
        int a = 1, b = 1;
        num = Math.abs(num);
        for (; a < limit && b < limit;) {
            if (a / (double) b > num) {
                b++;
            } else if (a / (double) b < num) {
                a++;
            } else {
                return b;
            }
        }
        return 0;
    }

    static public int rational_base(double num) {
        return rational_base(num, rational_limit);
    }

    static public boolean is_rational(double num, int limit) {
        return rational_base(num, limit) != 0;
    }

    static public boolean is_rational(double num) {
        return is_rational(num, rational_limit);
    }

    static public double max(double num1, double num2) {
        return (num1 > num2 ? num1 : num2);
    }

    static public int max(int num1, int num2) {
        return (num1 > num2 ? num1 : num2);
    }

    static public double min(double num1, double num2) {
        return (num1 > num2 ? num2 : num1);
    }

    static public int min(int num1, int num2) {
        return (num1 > num2 ? num2 : num1);
    }

    static public int gcf(int num1, int num2) {
        for (int i = min(num1, num2); i > 0; i--) {
            if (num1 % i == 0 && num2 % i == 0) {
                return i;
            }
        }
        return 1;
    }

    static public int gcf(int... nums) {
        if (nums.length < 2) {
            return 1;
        }
        int x = gcf(nums[0], nums[1]);
        for (int i = 2; i < nums.length; i++) {
            x = gcf(x, nums[i]);
        }
        return x;
    }

    static public int[] factor(int num) {
        boolean is_negative = num < 0;
        if (is_negative) {
            num *= 1;
        }
        int num_factors = 4; //we know any integer is divisible by 1 and itself, as well as -1 and its negative, so we know it has at least 4 factors
        int[] factors;
        if (num == 1 || num == -1) { //because 1 and 1 are the same number, they only count once towards the factors. same with -1
            factors = new int[2];
            factors[0] = 1;
            factors[1] = -1;
            return factors;
        }
        for (int i = 2; i < Math.sqrt(num); i++) {
            if (num % i == 0) {
                num_factors += 4; //the positive and negative multiples of the factors i and num/i
            }
        }
        factors = new int[num_factors];
        int j = 0;
        for (int i = 1; i < Math.sqrt(num); i++) {
            if (num % i == 0) {
                factors[j++] = i;
                factors[j++] = -1 * i;
                factors[j++] = num / i;
                factors[j++] = num / i * -1;
            }
        }
        return factors;
    }
}

class my_array_utils {

    static public Complex[] common_elements(Complex[] arr1, Complex[] arr2) {
        Complex[] comm_elem;
        int num_in_comm = 0;
        int k = 0;
        for (Complex arr11 : arr1) {
            for (Complex arr21 : arr2) {
                if (arr11.equals(arr21)) {
                    num_in_comm++;
                    break;
                }
            }
        }
        comm_elem = new Complex[num_in_comm];
        for (Complex arr11 : arr1) {
            for (Complex arr21 : arr2) {
                if (arr11.equals(arr21)) {
                    comm_elem[k++] = arr11;
                    break;
                }
            }
        }
        return comm_elem;
    }
}

class Complex {

    public double r;
    public double i;
    static public final Complex ZERO = new Complex(0, 0);

    Complex(double r, double i) {
        this.r = r;
        this.i = i;
    }

    Complex(double r) {
        this(r, 0);
    }

    Complex() {
        this(0, 0);
    }

    Complex(Complex c) {
        this(c.r, c.i);
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Complex)) {
            return false;
        }
        return (r == ((Complex) obj).r && i == ((Complex) obj).i);
    }

    @Override
    public String toString() {
        String str;
        if (i == 0) {
            if (my_math_utils.is_integer(r)) {
                str = String.format("%d", (long) r);
            } else {
                str = String.format("%f", r);
            }
        } else if (r == 0) {
            if (i == 1) {
                str = "i";
            } else {
                if (my_math_utils.is_integer(i)) {
                    str = String.format("%di", (long) i);
                } else {
                    str = String.format("%fi", i);
                }
            }
        } else {
            if (my_math_utils.is_integer(r)) {
                str = String.format("%d", (long) r);
            } else {
                str = String.format("%f", r);
            }
            if (my_math_utils.is_integer(i)) {
                str = str.concat(String.format("%di", (long) i));
            } else {
                str = str.concat(String.format("%fi", i));
            }
        }
        return str;
    }

    public Complex conj() {
        return new Complex(r, -1 * i);
    }

    public boolean isReal() {
        return i == 0.0;
    }

    public Complex add(Complex num) {
        return new Complex(num.r + r, num.i + i);
    }

    public Complex add(double num) {
        return new Complex(num + r, i);
    }

    public Complex multiply(Complex a) {
        return new Complex(r * a.r + -1 * i * a.i, r * a.i + i * a.r);
    }

    public Complex multiply(double a) {
        return new Complex(a * r, a * i);
    }

    static public Complex sqrt(double val) {
        if (val >= 0) {
            return new Complex(Math.sqrt(val));
        }
        return new Complex(0, Math.sqrt(-1 * val));
    }

    static public Complex sqrt(Complex val) {
        Complex root = new Complex();
        if (val.i == 0) {
            root = Complex.sqrt(val.r);
        } else {
            root.r = Math.sqrt((val.r + Math.sqrt(Math.pow(val.r, 2) + Math.pow(val.i, 2))) / 2.0);
            root.i = Math.signum(val.i) * Math.sqrt(((-1 * val.r) + Math.sqrt(Math.pow(val.r, 2) + Math.pow(val.i, 2))));
        }
        return root;
    }
}

class fraction {

    private boolean is_negative;
    public int numerator;
    public int denominator;
    public double value;

    fraction() {
        numerator = 0;
        denominator = 1;
        is_negative = false;
    }

    fraction(int val) {
        numerator = val;
        denominator = 1;
        is_negative = val < 0;
    }

    fraction(int num, int denom) {
        numerator = num;
        denominator = denom;
        is_negative = num < 0;
        if (denom < 0) {
            is_negative = !is_negative;
        }
    }

    fraction(double val) {
        value = val;
        numerator = 0;
        denominator = 1;
        is_negative = false;
        if (val < 0) {
            is_negative = true;
            val *= -1;
        }
        while (numerator / (double) denominator != val) {
            if (numerator / (double) denominator < val) {
                numerator++;
            } else {
                denominator++;
            }
        }
        if (is_negative) {
            numerator *= -1;
        }
    }

    @Override
    public String toString() {
        return String.format("(%d/%d)", numerator, denominator);
    }

    public void simplify() {
        for (int i = (numerator < denominator ? numerator : denominator); i > 0; i--) {
            if (numerator % i == 0 && denominator % i == 0) {
                numerator %= i;
                denominator %= i;
            }
        }
    }

    public double val() {
        return numerator / (double) denominator;
    }
}

class monomial {

    public int coefficient;
    public int exponent;

    monomial() {
        coefficient = 0;
        exponent = 0;
    }

    monomial(int val) {
        coefficient = val;
        exponent = 0;
    }

    monomial(int coeff, int expon) {
        coefficient = coeff;
        exponent = expon;
    }

    monomial(monomial mono) {
        this.coefficient = mono.coefficient;
        this.exponent = mono.exponent;
    }

    public polynomial add(monomial mono2) {
        polynomial sum;
        if (this.exponent == mono2.exponent) {
            sum = new polynomial(mono2);
            sum.elements[0].coefficient += this.coefficient;
        } else {
            sum = new polynomial(this, mono2);
        }
        return sum;
    }

    public monomial multiply(monomial mono2) {
        monomial product = new monomial();
        product.coefficient = this.coefficient * mono2.coefficient;
        product.exponent = this.exponent + mono2.exponent;
        return product;
    }

    public int compareto(monomial mono2) {
        if (this.exponent > mono2.exponent) {
            return 1;
        }
        if (this.exponent < mono2.exponent) {
            return -1;
        }
        return 0;
    }

    @Override
    public String toString() {
        String str;
        if (exponent == 0) {
            str = String.format("%d", coefficient);
        } else if (exponent == 1 && Math.abs(coefficient) > 1) {
            str = String.format("%dX", coefficient);
        } else if (exponent == 1 && coefficient == 1) {
            str = "X";
        } else if (exponent == 1 && coefficient == -1) {
            str = "-X";
        } else if (exponent > 1 && coefficient == 1) {
            str = String.format("X^%d", exponent);
        } else {
            str = String.format("%dX^%d", coefficient, exponent);
        }
        return str;
    }

    public double val(double x) {
        return (Math.pow(x, exponent) * coefficient);
    }
}

class polynomial {

    public int num_elements;
    public monomial[] elements;
    private boolean is_sorting;
    private boolean is_simplifying;

    polynomial() {
        num_elements = 1;
        elements = new monomial[1];
        elements[0] = new monomial();
    }

    polynomial(int num_elements) {
        this.num_elements = num_elements;
        elements = new monomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new monomial();
        }
    }

    polynomial(monomial... init_elements) {
        num_elements = init_elements.length;
        elements = new monomial[num_elements];
        int elements_initialized = 0;
        for (monomial temp_monomial : init_elements) {
            elements[elements_initialized++] = temp_monomial;
        }
    }

    polynomial(polynomial poly) {
        this.num_elements = poly.num_elements;
        this.elements = new monomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new monomial(poly.elements[i]);
        }
    }

    public int compareto(polynomial poly2) {
        if (this.num_elements < poly2.num_elements) {
            return 1;
        }
        if (this.num_elements > poly2.num_elements) {
            return -1;
        }
        if (this.elements[0].exponent > poly2.elements[0].exponent) {
            return 1;
        }
        if (this.elements[0].exponent < poly2.elements[0].exponent) {
            return -1;
        }
        return 0;
    }

    public boolean equals(polynomial poly2) {
        if (num_elements != poly2.num_elements) {
            return false;
        }
        for (int i = 0; i < num_elements; i++) {
            if (elements[i].coefficient != poly2.elements[i].coefficient || elements[i].exponent != poly2.elements[i].exponent) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String toString() {
        String str;
        if (num_elements == 1) {
            str = "";
        } else {
            str = "(";
        }
        str = str.concat(elements[0].toString());
        for (int i = 1; i < num_elements; i++) {
            if (elements[i].coefficient > 0) {
                str = str.concat(" +");
            }
            str = str.concat(" ".concat(elements[i].toString()));
        }
        if (num_elements != 1) {
            str = str.concat(")");
        }
        return str;
    }

    public void add_element() {
        monomial[] new_elements = new monomial[num_elements + 1];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        new_elements[num_elements++] = new monomial();
        elements = new_elements;
    }

    public void add_element(monomial new_element) {
        monomial[] new_elements = new monomial[num_elements + 1];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        new_elements[num_elements++] = new_element;
        elements = new_elements;
    }

    public void remove_element(int index) {
        monomial[] temp_elements = new monomial[num_elements - 1];
        System.arraycopy(this.elements, 0, temp_elements, 0, index);
        for (int i = index + 1; i < this.num_elements; i++) {
            temp_elements[i - 1] = this.elements[i];
        }
        this.elements = temp_elements;
        num_elements--;
    }

    public void simplify() {
        if (is_simplifying) {
            return;
        }
        if (num_elements == 1) {
            return;
        }
        is_simplifying = true;
        for (int i = 0; i < num_elements; i++) {
            for (int j = i + 1; j < num_elements;) {
                if (elements[i].exponent == elements[j].exponent) {
                    elements[i].coefficient += elements[j].coefficient;
                    this.remove_element(j);
                } else {
                    j++;
                }
            }
        }
        this.sort();
        is_simplifying = false;
    }

    public void sort() {
        if (is_sorting) {
            return;
        }
        is_sorting = true;
        for (int i = this.num_elements - 1; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                if (this.elements[j].exponent < this.elements[j + 1].exponent) {
                    monomial temp_mono = this.elements[j];
                    this.elements[j] = this.elements[j + 1];
                    this.elements[j + 1] = temp_mono;
                }
            }
        }
        this.simplify();
        is_sorting = false;
    }

    public polynomial add(polynomial poly2) {
        polynomial sum = new polynomial(this.num_elements + poly2.num_elements);
        for (int i = 0; i < this.num_elements; i++) {
            sum.elements[i] = new monomial(this.elements[i]);
        }
        for (int i = this.num_elements; i < this.num_elements + poly2.num_elements; i++) {
            sum.elements[i] = new monomial(poly2.elements[i - this.num_elements]);
        }
        this.simplify();
        return sum;
    }

    public polynomial multiply(polynomial poly2) {
        polynomial product = new polynomial(this.num_elements * poly2.num_elements);
        for (int i = 0; i < product.num_elements; i++) {
            product.elements[i] = new monomial(this.elements[i % this.num_elements].multiply(poly2.elements[i / this.num_elements]));
        }
        product.simplify();
        return product;
    }

    public Complex[] find_zeros() {
        this.sort();
        Complex[] zeros = new Complex[elements[0].exponent];
        for (int i = 0; i < zeros.length; i++) {
            zeros[i] = new Complex();
        }
        int num_zeros = 0;
        switch (num_elements) {
            case 1:
                for (Complex zero : zeros) {
                    zero.r = 0.0;
                }
                return zeros;
            case 2:
                if (this.elements[0].exponent == 1) {
                    zeros[0].r = ((-1 * elements[1].coefficient) / (double) elements[0].coefficient);
                    return zeros;
                } else if (elements[0].exponent % 2 == 0) {
                    zeros[0] = Complex.sqrt(-1 * elements[1].coefficient / (double) elements[0].coefficient);
                    zeros[1] = zeros[0].conj();
                    return zeros;
                } else if (elements[0].exponent % 3 == 0) {
                    zeros[0].r = zeros[1].r = zeros[2].r = Math.cbrt(elements[1].coefficient / (double) elements[0].coefficient);
                    return zeros;
                }
            case 3:
                /*
                 *  use quadratic formula
                 *      -b +- sqrt(b^2-4ac)
                 *  x = __________________
                 *              2a
                 */
                double a = elements[0].coefficient,
                 b = elements[1].coefficient,
                 c = elements[2].coefficient;
                zeros[0] = Complex.sqrt(b * b - 4 * a * c).add(-1 * b).multiply(0.5 / a);
                zeros[1] = Complex.sqrt(b * b - 4 * a * c).multiply(-1).add(-1 * b).multiply(0.5 / a);
                break;
            default:
                /*
                 add methods recently covered in precalc to find zeros of longer polynomials
                 */
                double[] possible_roots;
                int[] trailing_factors = my_math_utils.factor(elements[num_elements - 1].coefficient);
                int[] leading_factors = my_math_utils.factor(elements[0].coefficient);
                possible_roots = new double[trailing_factors.length * leading_factors.length];
                int k = 0;
                for (int i = 0; i < trailing_factors.length; i++) {
                    for (int j = 0; j < leading_factors.length; j++) {
                        possible_roots[k++] = trailing_factors[i] / (double) leading_factors[j];
                    }
                }
                k = 0;
                for (int i = 0; i < possible_roots.length; i++) {
                    if (this.val(possible_roots[i]) == 0.0) {
                        zeros[k++].r = possible_roots[i];
                    }
                }
                java.util.Arrays.sort(trailing_factors);
                for (int i = -100; i < 100; i++) {
                    if (java.util.Arrays.binarySearch(trailing_factors, i) != 0) {
                        continue;
                    }
                    if (val(i) == 0.0) {
                        zeros[k++].r = val(i);
                    }
                }
        }
        for (; num_zeros < zeros.length; num_zeros++) {
            if (zeros[num_zeros].r == 0.0 || Double.isNaN(zeros[num_zeros].r) || Double.isInfinite(zeros[num_zeros].r)) {
                num_zeros = my_math_utils.max(num_zeros - 1, 0);
                break;
            }
        }
        Complex[] real_zeros = new Complex[num_zeros];
        System.arraycopy(zeros, 0, real_zeros, 0, num_zeros);
        return real_zeros;
    }

    public double val(double x) {
        double sum = 0;
        for (int i = 0; i < num_elements; i++) {
            sum += elements[i].val(x);
        }
        return sum;
    }

    public int can_factor() {
        int[] factors1;
        int[] factors2;
        double[] possible_roots;
        int k;
        this.sort();
        if (num_elements > 1 && elements[num_elements - 1].exponent != 0) {
            return 1;
        }
        if (elements[0].coefficient < 0) {
            return -1;
        }
        switch (this.num_elements) {
            case 1:
                return 0;
            case 2:
                if (elements[0].exponent == 1) {
                    return 0;
                }
                if (elements[0].exponent % 2 == 0) {
                    if (my_math_utils.is_square(this.elements[0].coefficient)) {
                        if (my_math_utils.is_negative(this.elements[1].coefficient)) {
                            if (my_math_utils.is_square((-1 * this.elements[1].coefficient))) {
                                return 2;
                            }
                        }
                    }
                }
                if (elements[0].exponent % 3 == 0) {
                    if (my_math_utils.is_cube(this.elements[0].coefficient)) {
                        if (my_math_utils.is_cube(this.elements[1].coefficient)) {
                            return 3;
                        }
                    }
                }
                break;
            case 3:
                if (this.elements[0].exponent - this.elements[1].exponent != this.elements[1].exponent - this.elements[2].exponent) {
                    return 0;
                }
                if (my_math_utils.is_negative(this.elements[0].coefficient)) {
                    for (int i = this.elements[0].coefficient; i <= -1 * this.elements[0].coefficient; i--) {
                        if (i == 0) {
                            continue;
                        }
                        if (this.elements[0].coefficient % i != 0) {
                            continue;
                        }
                        if (my_math_utils.is_negative(this.elements[2].coefficient)) {
                            for (int j = this.elements[2].coefficient; j <= -1 * this.elements[2].coefficient; j++) {
                                if (j == 0) {
                                    continue;
                                }
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        } else {
                            for (int j = this.elements[2].coefficient; j >= 1; j--) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                            for (int j = -1 * this.elements[2].coefficient; j <= -1; j++) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        }
                    }
                } else {
                    for (int i = 1; i <= this.elements[0].coefficient; i++) {
                        if (this.elements[0].coefficient % i != 0) {
                            continue;
                        }
                        if (my_math_utils.is_negative(this.elements[2].coefficient)) {
                            for (int j = this.elements[2].coefficient; j <= -1 * this.elements[2].coefficient; j++) {
                                if (j == 0) {
                                    continue;
                                }
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        } else {
                            for (int j = this.elements[2].coefficient; j >= 1; j--) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                            for (int j = -1 * this.elements[2].coefficient; j <= -1; j++) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        }
                    }
                    for (int i = -1; i <= -1 * this.elements[0].coefficient; i--) {
                        if (this.elements[0].coefficient % i != 0) {
                            continue;
                        }
                        if (my_math_utils.is_negative(this.elements[2].coefficient)) {
                            for (int j = this.elements[2].coefficient; j <= -1 * this.elements[2].coefficient; j++) {
                                if (j == 0) {
                                    continue;
                                }
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        } else {
                            for (int j = this.elements[2].coefficient; j >= 1; j--) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                            for (int j = -1 * this.elements[2].coefficient; j <= -1; j++) {
                                if (this.elements[2].coefficient % j != 0) {
                                    continue;
                                }
                                if (j * (this.elements[0].coefficient / i) + i * (this.elements[2].coefficient / j) == this.elements[1].coefficient) {
                                    return (j < 0 ? 100000 : 0) + Math.abs(j) * 1000 + (i < 0 ? 100 : 0) + Math.abs(i);
                                }
                            }
                        }
                    }
                }
                return 0;
            case 4:
                if (elements[0].coefficient / (double) elements[2].coefficient == elements[1].coefficient / (double) elements[3].coefficient) {
                    if (elements[0].coefficient / (double) elements[1].coefficient == elements[2].coefficient / (double) elements[3].coefficient) {
                        if (elements[0].exponent - elements[2].exponent == elements[1].exponent - elements[3].exponent) {
                            if (elements[0].exponent - elements[1].exponent == elements[2].exponent - elements[3].exponent) {
                                return 4;
                            }
                        }
                    }
                }
                /*
                 test if it can be broken into a binomial and a trinomial by using the rational root theorem and testing the results using synthetic substitution
                 */
                factors1 = my_math_utils.factor(elements[0].coefficient);
                factors2 = my_math_utils.factor(elements[3].coefficient);
                possible_roots = new double[factors1.length * factors2.length / 2];
                k = 0;
                for (int i = 0; i < factors1.length; i += 2) {
                    for (int j = 0; j < factors2.length; j++) {
                        possible_roots[k++] = factors1[i] / (double) factors2[j];
                    }
                }
                for (int i = 0; i < possible_roots.length; i++) {
                    if (val(possible_roots[i]) == 0) {
                        int a, b;
                        a = my_math_utils.rational_base(possible_roots[i]);
                        b = (int) (possible_roots[i] * a);
                        return a * 10000000 + b;
                    }
                }
                break;
            default:
                //use the rational root theorem to find the possible real roots, then test those roots using synthetic substitution
                factors1 = my_math_utils.factor(elements[0].coefficient);
                factors2 = my_math_utils.factor(elements[num_elements - 1].coefficient);
                possible_roots = new double[factors1.length * factors2.length];
                k = 0;
                for (int i = 0; i < factors1.length; i++) {
                    for (int j = 0; j < factors2.length; j++) {
                        possible_roots[k++] = factors1[i] / (double) factors2[j];
                    }
                }
                for (int i = 0; i < possible_roots.length; i++) {
                    if (val(possible_roots[i]) == 0.0) {
                        int a, b;
                        a = my_math_utils.rational_base(possible_roots[i], Integer.MAX_VALUE / 10000000 - 1); //not sure if subtracting one is neccessary, but I'm doing it just to be sure :)
                        b = (int) (possible_roots[i] * a);
                        return a * 10000000 + b;
                    }
                }
                double v;
                for (int i = -100; i < 100; i++) {
                    v = val(i);
                    System.out.printf("f(%d) = %f\n", i, v);
                    if (v == 0.0) {
                        return 10000000 + i;
                    }
                }
        }
        return 0;
    }

    /*
     * returns true if argument is a factor, otherwise returns false.
     * if argument is a factor, then the current instance will be divided by the argument
     * if argument is a linear binomial, it is advised to test it first with val() before trying to divide by it; in very large expressions, this may be more resource-efficient
     */
    public boolean long_division(polynomial divisor) {
        boolean is_successful;
        if (divisor.num_elements < divisor.elements[0].exponent + 1) {
            polynomial expanded_divisor = new polynomial(divisor.elements[0].exponent + 1);
            for (int i = 0; i < expanded_divisor.num_elements; i++) {
                expanded_divisor.elements[i].exponent = expanded_divisor.num_elements - i - 1;
            }
            for (int i = 0; i < divisor.num_elements; i++) {
                expanded_divisor.elements[expanded_divisor.num_elements - divisor.elements[i].exponent - 1].coefficient = divisor.elements[i].coefficient;
            }
            divisor = expanded_divisor;
        }
        polynomial expanded_elements;
        if (num_elements < elements[0].exponent) {
            expanded_elements = new polynomial(elements[0].exponent);
            for (int i = 0; i < expanded_elements.num_elements; i++) {
                expanded_elements.elements[i].exponent = expanded_elements.num_elements - i - 1;
            }
            for (int i = 0; i < num_elements; i++) {
                expanded_elements.elements[expanded_elements.num_elements - elements[i].exponent /* - 1*/].coefficient = elements[i].coefficient;
            }
        } else {
            expanded_elements = new polynomial(this);
        }
        double[] dividend = new double[expanded_elements.num_elements];
        for (int i = 0; i < dividend.length; i++) {
            dividend[i] = expanded_elements.elements[i].coefficient;
        }
        double[] quotient = new double[expanded_elements.num_elements];
        for (int i = 0; i < quotient.length - divisor.num_elements + 1; i++) {
            quotient[i] = dividend[i] / divisor.elements[0].coefficient;
            for (int j = 0; j < divisor.num_elements; j++) {
                dividend[i + j] -= quotient[i] * divisor.elements[j].coefficient;
            }
        }
        is_successful = quotient[quotient.length - 1] == 0.0;
        if (is_successful) {
            for (int i = 0; i < quotient.length; i++) {
                if (!my_math_utils.is_integer(quotient[i])) { //for the sake of simplicity, after it divides, if any coefficient is not an integer, it will treat it as now being divisible, even if the remainder is 0
                    return false;
                }
            }
            if (expanded_elements.num_elements == num_elements) {
                this.remove_element(0);
                for (int i = 0; i < num_elements; i++) {
                    elements[i].coefficient = (int) quotient[i];
                }
            } else {
                expanded_elements.remove_element(0);
                for (int i = 0; i < num_elements; i++) {
                    expanded_elements.elements[i].coefficient = (int) quotient[i];
                }
                elements = expanded_elements.elements;
                num_elements = expanded_elements.num_elements;
            }
        }
        return is_successful;
    }

    static public polynomial add(polynomial poly1, polynomial poly2) {
        return poly1.add(poly2);
    }

    static public polynomial multiply(polynomial poly1, polynomial poly2) {
        return poly1.multiply(poly2);
    }
}

class poly_n_d {

    public int num_elements;
    public polynomial[] elements;
    private boolean is_sorting;
    private boolean is_simplifying;

    poly_n_d() {
        num_elements = 1;
        elements = new polynomial[1];
        elements[0] = new polynomial();
    }

    poly_n_d(int num_elements) {
        this.num_elements = num_elements;
        elements = new polynomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            this.elements[i] = new polynomial();
        }
    }

    poly_n_d(polynomial... init_elements) {
        this.num_elements = init_elements.length;
        elements = new polynomial[this.num_elements];
        int elements_initialized = 0;
        for (polynomial temp_polynomial : init_elements) {
            elements[elements_initialized++] = temp_polynomial;
        }
    }

    poly_n_d(poly_n_d poly) {
        num_elements = poly.num_elements;
        elements = new polynomial[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new polynomial(poly.elements[i]);
        }
    }

    @Override
    public String toString() {
        this.sort();
        String str = elements[0].toString();
        for (int i = 1; i < num_elements; i++) {
            str = str.concat(" ".concat(elements[i].toString()));
        }
        return str;
    }

    public double val(double x) {
        double sum = 1;
        for (int i = 0; i < num_elements; i++) {
            sum *= elements[i].val(x);
        }
        return sum;
    }

    public void add_element() {
        polynomial[] new_elements = new polynomial[num_elements + 1];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        new_elements[num_elements] = new polynomial();
        num_elements += 1;
        elements = new_elements;
    }

    public void add_element(polynomial new_element) {
        polynomial[] new_elements = new polynomial[num_elements + 1];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        new_elements[num_elements] = new_element;
        elements = new_elements;
        num_elements += 1;
    }

    public void add_element(int num_to_add) {
        polynomial[] new_elements = new polynomial[num_elements + num_to_add];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        for (int i = 0; i < num_to_add; i++) {
            new_elements[num_elements + i] = new polynomial();
        }
        num_elements += num_to_add;
        elements = new_elements;
    }

    public void remove_element(int index) {
        polynomial[] temp_elements = new polynomial[num_elements - 1];
        System.arraycopy(this.elements, 0, temp_elements, 0, index);
        for (int i = index + 1; i < this.num_elements; i++) {
            temp_elements[i - 1] = this.elements[i];
        }
        this.elements = temp_elements;
        num_elements--;
    }

    public void sort() {
        if (!is_simplifying) {
            simplify();
        }
        if (is_sorting) {
            return;
        }
        is_sorting = true;
        for (int i = 0; i < num_elements; i++) {
            elements[i].sort();
        }
        for (int i = num_elements - 1; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                if (elements[j].num_elements < elements[j + 1].num_elements) {
                    polynomial temp_poly = elements[j];
                    elements[j] = elements[j + 1];
                    elements[j + 1] = temp_poly;
                }
            }
        }
        simplify();
        is_sorting = false;
    }

    public void simplify() {
        if (is_simplifying) {
            return;
        }
        if (num_elements == 1) {
            return;
        }
        is_simplifying = true;
        for (int i = 0; i < num_elements; i++) {
            if (elements[i].num_elements == 1) {
                for (int j = i + 1; j < num_elements;) {
                    if (elements[j].num_elements == 1) {
                        elements[i] = elements[i].multiply(elements[j]);
                        this.remove_element(j);
                        continue;
                    }
                    j++;
                }
            }
        }
        sort();
        is_simplifying = false;
    }

    public void factor() {
        this.sort();
        int lowest_exponent;
        int factorability;
        for (int i = 0; i < this.num_elements;) {
            if ((factorability = elements[i].can_factor()) == 0) {
                i++;
                continue;
            }
            if (factorability == -1) {
                this.add_element(new polynomial(new monomial(-1, 0)));
                for (int j = 0; j < elements[i].num_elements; j++) {
                    elements[i].elements[j].coefficient *= -1;
                }
            } else if (factorability == 1) {
                int poly_gcf = elements[i].elements[0].coefficient;
                for (int j = 1; j < elements[i].num_elements; j++) {
                    poly_gcf = my_math_utils.gcf(poly_gcf, elements[i].elements[j].coefficient);
                }
                this.add_element(new polynomial(new monomial(poly_gcf, (lowest_exponent = this.elements[i].elements[this.elements[i].num_elements - 1].exponent))));
                for (int j = 0; j < elements[i].num_elements; j++) {
                    elements[i].elements[j].coefficient /= poly_gcf;
                    elements[i].elements[j].exponent -= lowest_exponent;
                }
            } else if (factorability == 2) {
                elements[i].elements[0].coefficient = (int) Math.round(Math.sqrt(elements[i].elements[0].coefficient));
                elements[i].elements[0].exponent /= 2;
                elements[i].elements[1].coefficient = (int) Math.round(Math.sqrt(-1 * elements[i].elements[1].coefficient));
                elements[i].elements[1].exponent /= 2;
                this.add_element(new polynomial(new monomial(elements[i].elements[0]), new monomial(elements[i].elements[1])));
                elements[num_elements - 1].elements[1].coefficient *= -1;
            } else if (factorability == 3) {
                if (elements[i].elements[1].coefficient < 0) {
                    elements[i].elements[0].coefficient = (int) Math.cbrt(elements[i].elements[0].coefficient);
                    elements[i].elements[0].exponent /= 3;
                    elements[i].elements[1].coefficient = (int) Math.cbrt(elements[i].elements[1].coefficient);
                    elements[i].elements[1].exponent /= 3;
                    this.add_element(new polynomial(3));
                    int j = num_elements - 1;
                    elements[j].elements[0].coefficient = (int) Math.pow(elements[i].elements[0].coefficient, 2);
                    elements[j].elements[0].exponent = 2 * elements[i].elements[0].exponent;
                    elements[j].elements[1].coefficient = -1 * elements[i].elements[0].coefficient * elements[i].elements[1].coefficient;
                    elements[j].elements[1].exponent = elements[i].elements[0].exponent;
                    elements[j].elements[2].coefficient = (int) Math.pow(elements[i].elements[1].coefficient, 2);
                } else {
                    elements[i].elements[0].coefficient = (int) Math.cbrt(elements[i].elements[0].coefficient);
                    elements[i].elements[0].exponent /= 3;
                    elements[i].elements[1].coefficient = (int) Math.cbrt(elements[i].elements[1].coefficient);
                    elements[i].elements[1].exponent /= 3;
                    this.add_element(new polynomial(3));
                    int j = num_elements - 1;
                    elements[j].elements[0].coefficient = (int) Math.pow(elements[i].elements[0].coefficient, 2);
                    elements[j].elements[0].exponent = 2 * elements[i].elements[0].exponent;
                    elements[j].elements[1].coefficient = -1 * elements[i].elements[0].coefficient * elements[i].elements[1].coefficient;
                    elements[j].elements[1].exponent = elements[i].elements[0].exponent;
                    elements[j].elements[2].coefficient = -1 * ((int) Math.pow(elements[i].elements[1].coefficient, 2));
                }
            } else if (factorability == 4) { //x^3+6x^2+4x+24   ->      (x^2+4)(x+6)        (aX^2+b)(cX+d) = acX^3 + adX^2 + bcX + bd
                // http://www.wikihow.com/Factor-a-Cubic-Polynomial
                this.add_element(new polynomial(2));
                int j = num_elements - 1;
                elements[j].elements[0].coefficient = my_math_utils.gcf(elements[i].elements[0].coefficient, elements[i].elements[1].coefficient);
                elements[j].elements[0].exponent = elements[i].elements[1].exponent;
                elements[j].elements[1].coefficient = my_math_utils.gcf(Math.abs(elements[i].elements[2].coefficient), Math.abs(elements[i].elements[3].coefficient))
                        * (elements[i].elements[3].coefficient < 0 ? -1 : 1);
                elements[j].elements[1].exponent = 0;
                this.add_element(new polynomial(2));
                j++;
                elements[j].elements[0].coefficient = elements[i].elements[0].coefficient / elements[j - 1].elements[0].coefficient;
                elements[j].elements[0].exponent = elements[i].elements[0].exponent - elements[i].elements[1].exponent;
                elements[j].elements[1].coefficient = elements[i].elements[1].coefficient / elements[j - 1].elements[0].coefficient;
                elements[j].elements[1].exponent = 0;
                this.remove_element(i);
            } else if (factorability >= 1000 && factorability < 9000000) {
                int a = (factorability % 100) * ((factorability / 100) % 10 == 1 ? -1 : 1);
                int b = ((factorability / 1000) % 100) * (factorability / 100000 == 1 ? -1 : 1);
                this.add_element(new polynomial(2));
                int j = num_elements - 1;
                elements[j].elements[0].coefficient = a;
                elements[j].elements[0].exponent = elements[i].elements[0].exponent / 2;
                elements[j].elements[1].coefficient = b;
                elements[j].elements[1].exponent = 0;
                this.add_element(new polynomial(2));
                j++;
                elements[j].elements[0].coefficient = elements[i].elements[0].coefficient / a;
                elements[j].elements[0].exponent = elements[i].elements[0].exponent / 2;
                elements[j].elements[1].coefficient = elements[i].elements[2].coefficient / b;
                elements[j].elements[1].exponent = 0;
                this.remove_element(i);
            } else if (Math.abs(factorability) >= 9000000) { //manual factoring out of a linear binomial found to be a factor in can_factor()
                int a = Math.round(factorability / (float) 10000000);
                int b = factorability - (a * 10000000);
                polynomial factor = new polynomial(2);
                factor.elements[0].coefficient = a;
                factor.elements[0].exponent = 1;
                factor.elements[1].coefficient = -1 * b;
                if (!elements[i].long_division(factor)) {
                    System.out.println("Error factoring: results of long division and synthetic substitution do not match");
                    System.exit(1);
                }
                this.add_element(factor);
            } else {

            }
            this.sort();
            i = 0;
        }
    }

    public Complex[] find_zeros() {
        int num_zeros = 0;
        int curr_zero = 0;
        Complex[] zeros;
        Complex[][] element_zeros = new Complex[num_elements][];
        for (int i = 0; i < num_elements; i++) {
            element_zeros[i] = elements[i].find_zeros();
        }
        for (Complex[] element_zero : element_zeros) {
            num_zeros += element_zero.length;
        }
        zeros = new Complex[num_zeros];
        for (int i = 0; i < num_elements; i++) {
            for (Complex element_zero : element_zeros[i]) {
                zeros[curr_zero++] = element_zero;
            }
        }
        return zeros;
    }

    public void distribute() {
        while (num_elements > 1) {
            elements[0] = elements[0].multiply(elements[1]);
            this.remove_element(1);
        }
    }

    static public poly_n_d multiply(poly_n_d poly1, poly_n_d poly2) {
        poly_n_d product = new poly_n_d(poly1);
        for (int i = 0; i < poly2.num_elements; i++) {
            product.add_element(poly2.elements[i]);
        }
        return product;
    }

    static public poly_n_d add(poly_n_d poly1, poly_n_d poly2) {
        poly_n_d sum = new poly_n_d(poly1);
        sum.distribute();
        poly_n_d temp_poly = new poly_n_d(poly2);
        temp_poly.distribute();
        sum.elements[0] = sum.elements[0].add(temp_poly.elements[0]);
        return sum;
    }
}

class polynomial_fraction {

    public poly_n_d numerator;
    public poly_n_d denominator;

    polynomial_fraction() {
        numerator = new poly_n_d();
        denominator = new poly_n_d();
        denominator.elements[0].elements[0].coefficient = 1;
    }

    polynomial_fraction(poly_n_d numerator) {
        this.numerator = numerator;
        this.denominator = new poly_n_d();
        this.denominator.elements[0].elements[0].coefficient = 1;
    }

    polynomial_fraction(poly_n_d numerator, poly_n_d denominator) {
        this.numerator = numerator;
        this.denominator = denominator;
    }

    polynomial_fraction(polynomial_fraction poly_fract) {
        numerator = new poly_n_d(poly_fract.numerator);
        denominator = new poly_n_d(poly_fract.denominator);
    }

    public double val(double x) {
        return numerator.val(x) / denominator.val(x);
    }

    public void factor() {
        numerator.factor();
        denominator.factor();
    }

    public Complex[][] find_solutions() {
        Complex[][] solutions = new Complex[5][];
        Complex[] zeros;
        Complex[] holes;
        Complex[] h_asymptotes;
        Complex[] v_asymptotes;
        int num_degree = 0;
        for (int i = 0; i < numerator.num_elements; i++) {
            num_degree += numerator.elements[i].elements[0].exponent;
        }
        int denom_degree = 0;
        for (int i = 0; i < denominator.num_elements; i++) {
            denom_degree += denominator.elements[i].elements[0].exponent;
        }
        if (num_degree > denom_degree) {
            h_asymptotes = new Complex[0];
        } else if (num_degree < denom_degree) {
            h_asymptotes = new Complex[1];
            h_asymptotes[0] = new Complex();
        } else {
            h_asymptotes = new Complex[1];
            h_asymptotes[0] = new Complex();
            int num_coeff = 1;
            for (int i = 0; i < numerator.num_elements; i++) {
                num_coeff *= numerator.elements[i].elements[0].coefficient;
            }
            int denom_coeff = 1;
            for (int i = 0; i < denominator.num_elements; i++) {
                denom_coeff *= denominator.elements[i].elements[0].coefficient;
            }
            h_asymptotes[0].r = num_coeff / (double) denom_coeff;
        }
        Complex[] num_zeros, denom_zeros;
        num_zeros = numerator.find_zeros();
        denom_zeros = denominator.find_zeros();
        holes = my_array_utils.common_elements(num_zeros, denom_zeros);
        zeros = new Complex[num_zeros.length - holes.length];
        for (int i = 0; i < zeros.length; i++) {
            zeros[i] = new Complex();
        }
        v_asymptotes = new Complex[denom_zeros.length - holes.length];
        for (int i = 0; i < v_asymptotes.length; i++) {
            v_asymptotes[i] = new Complex();
        }
        int j = 0;
        int k = 0;
        if (holes.length > 0) {
            for (int i = 0; i < num_zeros.length; i++) {
                if (num_zeros[i].equals(holes[j])) {
                    j++;
                } else {
                    zeros[k++] = num_zeros[i];
                }
                if (j == holes.length) {
                    for (; i < num_zeros.length && k < zeros.length; i++) {
                        zeros[k++] = num_zeros[i];
                    }
                }
            }
            j = 0;
            k = 0;
            for (int i = 0; i < denom_zeros.length; i++) {
                if (denom_zeros[i].equals(holes[j])) {
                    j++;
                } else {
                    v_asymptotes[k++] = denom_zeros[i];
                }
                if (j == holes.length) {
                    for (; i < denom_zeros.length && k < v_asymptotes.length; i++) {
                        v_asymptotes[k++] = denom_zeros[i];
                    }
                }
            }
            solutions[0] = zeros;
            solutions[3] = v_asymptotes;
        } else {
            solutions[0] = num_zeros;
            solutions[3] = denom_zeros;
        }
        solutions[1] = holes;
        solutions[2] = h_asymptotes;
        return solutions;
    }

    @Override
    public String toString() {
        String str = numerator.toString();
        if (!(denominator.num_elements == 1 && denominator.elements[0].num_elements == 1 && denominator.elements[0].elements[0].coefficient == 1)) {
            str = str.concat(" / ");
            str = str.concat(denominator.toString());
        }

        /* *
         String num_str = numerator.toString();
         int num_length = num_str.length();
         String denom_str = denominator.toString();
         int denom_length = denom_str.length();
         if (num_length > denom_length) {
         str = num_str.concat("\n");
         for (int i = 0; i < num_length; i++) {
         str = str.concat("-");
         }
         str = str.concat("\n");
         for (int i = 0; i < Math.floor((num_length - denom_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         str = str.concat(denom_str);
         for (int i = 0; i < Math.ceil((num_length - denom_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         } else {
         for (int i = 0; i < Math.floor((denom_length - num_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         str = str.concat(num_str);
         for (int i = 0; i < Math.ceil((denom_length - num_length) / (double) 2); i++) {
         str = str.concat(" ");
         }
         str = num_str.concat("\n");
         for (int i = 0; i < denom_length; i++) {
         str = str.concat("-");
         }
         str = str.concat("\n");
         str = str.concat(denom_str);
         }
         /* */
        return str;
    }

    public polynomial_fraction add(polynomial_fraction poly_fract2) {
        polynomial_fraction sum = new polynomial_fraction();
        sum.denominator = poly_n_d.multiply(denominator, poly_fract2.denominator);
        sum.numerator = poly_n_d.add(poly_n_d.multiply(numerator, poly_fract2.denominator), poly_n_d.multiply(poly_fract2.numerator, denominator));
        return sum;
    }
}

class expression {

    public int num_elements;
    public polynomial_fraction[] elements;

    expression() {
        num_elements = 1;
        elements = new polynomial_fraction[1];
        elements[0] = new polynomial_fraction();
    }

    expression(int num_elements) {
        this.num_elements = num_elements;
        elements = new polynomial_fraction[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new polynomial_fraction();
        }
    }

    expression(polynomial_fraction... init_elements) {
        this.num_elements = init_elements.length;
        elements = new polynomial_fraction[num_elements];
        int elements_initialized = 0;
        for (polynomial_fraction temp_poly_fract : init_elements) {
            elements[elements_initialized++] = temp_poly_fract;
        }
    }

    expression(expression expres) {
        this.num_elements = expres.num_elements;
        this.elements = new polynomial_fraction[num_elements];
        for (int i = 0; i < num_elements; i++) {
            elements[i] = new polynomial_fraction(expres.elements[i]);
        }
    }

    public void add_element() {
        polynomial_fraction[] new_elements = new polynomial_fraction[num_elements + 1];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        new_elements[num_elements] = new polynomial_fraction();
        elements = new_elements;
    }

    public void add_element(polynomial_fraction new_element) {
        polynomial_fraction[] new_elements = new polynomial_fraction[num_elements + 1];
        System.arraycopy(elements, 0, new_elements, 0, num_elements);
        new_elements[num_elements] = new_element;
        elements = new_elements;
    }

    public void remove_element(int element_index) {
        polynomial_fraction[] new_elements = new polynomial_fraction[num_elements - 1];
        System.arraycopy(elements, 0, new_elements, 0, element_index);
        for (int i = element_index + 1; i < this.num_elements; i++) {
            new_elements[i - 1] = elements[i];
        }
        elements = new_elements;
    }

    public double val(double x) {
        double sum = 0;
        for (int i = 0; i < num_elements; i++) {
            sum += elements[i].val(x);
        }
        return sum;
    }

    public void factor() {
        for (int i = 0; i < num_elements; i++) {
            elements[i].factor();
        }
    }

    @Override
    public String toString() {
        String str = elements[0].toString();
        for (int i = 1; i < num_elements; i++) {
            str = str.concat(elements[i].toString());
        }
        return str;
    }
}

class equation {

    public expression left_side;
    public expression right_side;

    equation() {
        left_side = new expression();
        right_side = new expression();
    }

    equation(expression left_side) {
        this.left_side = left_side;
        this.right_side = new expression();
    }

    equation(expression left_side, expression right_side) {
        this.left_side = left_side;
        this.right_side = right_side;
    }

    equation(String input_str) {
        this.left_side = new expression();
        this.right_side = new expression();
        parse_input(input_str);
    }

    public final boolean parse_input(String input_str) {
        boolean has_equals = false;
        int curr_monomial = 0;
        int curr_polynomial = 0;
        int curr_poly_fract = 0;
        final int TOP = 0;
        final int BOTTOM = 1;
        int curr_fract_part = TOP;
        final int MONOMIAL = 0;
        final int POLYNOMIAL = 1;
        final int POLY_N_D = 2;
        final int POLY_FRACT = 3;
        final int EXPRESSION = 4;
        int curr_state = EXPRESSION;
        boolean is_negative = false;
        final int COEFFICIENT = 0;
        final int EXPONENT = 1;
        int curr_num = COEFFICIENT;
        int ch;
        char variable_char = 0;
        int length = input_str.length();
        for (int i = 0; i < length; i++) {
            if (input_str.charAt(i) == '=') {
                has_equals = true;
            }
        }
        for (int i = 0; i < length && input_str.charAt(i) != '='; i++) {
            ch = input_str.charAt(i);
            if (ch == '(') {
                if (curr_num == EXPONENT) {
                    System.out.println("error: cannot have non-constant exponent");
                    return false;
                }
                curr_monomial = 0;
                if (curr_state == EXPRESSION || curr_state == POLY_FRACT) {
                    for (int j = i + 1; j < length; j++) {
                        if (input_str.charAt(j) == ')') {
                            break;
                        }
                        if (j == length - 1) {
                            System.out.println("error: unpaired open parenthesis");
                            return false;
                        }
                    }
                } else {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.add_element();
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.add_element();
                    }
                    curr_polynomial++;
                }
                curr_state = POLYNOMIAL;
            } else if (ch == ')') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                    is_negative = false;
                }
                if (curr_state == POLY_N_D) {
                    continue;
                }
                curr_state = POLY_N_D;
                curr_monomial = 0;
                curr_num = COEFFICIENT;
            } else if (ch == '+') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                    is_negative = false;
                }
                if (curr_state == MONOMIAL) {
                    curr_num = COEFFICIENT;
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                    }
                    curr_monomial++;
                } else if (curr_state == POLY_N_D) {
                    this.left_side.add_element();
                    curr_poly_fract++;
                    curr_num = COEFFICIENT;
                    curr_polynomial = 0;
                    curr_fract_part = TOP;
                    curr_monomial = 0;
                }
            } else if (ch == '-') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                    is_negative = false;
                }
                if (curr_state == MONOMIAL) {
                    curr_num = COEFFICIENT;
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                    }
                    curr_monomial++;
                    is_negative = true;
                } else if (curr_state == POLYNOMIAL) {
                    is_negative = true;
                } else if (curr_state == POLY_N_D) {
                    this.left_side.add_element();
                    curr_poly_fract++;
                    curr_num = COEFFICIENT;
                    curr_polynomial = 0;
                    curr_fract_part = TOP;
                    curr_monomial = 0;
                    is_negative = true;
                }
            } else if (ch == '/') {
                if (is_negative) {
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    }
                }
                if (curr_state == MONOMIAL || curr_state == POLY_N_D) {
                    curr_state = POLY_FRACT;
                }
                curr_monomial = 0;
                curr_polynomial = 0;
                curr_fract_part = BOTTOM;
                curr_num = COEFFICIENT;
                this.left_side.elements[curr_poly_fract].denominator.elements[0].elements[0].coefficient = 0;
                is_negative = false;
            } else if (ch >= '0' && ch <= '9') {
                curr_state = MONOMIAL;
                if (curr_fract_part == TOP) {
                    if (curr_num == COEFFICIENT) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                    } else {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                    }
                } else {
                    if (curr_num == COEFFICIENT) {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                    }
                }
            } else if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')) {
                if (curr_num == EXPONENT) {
                    System.out.println("error: cannot have non-constant exponent");
                    return false;
                }
                if (ch >= 'a' && ch <= 'z') {
                    ch += 'A' - 'a';
                }
                if (variable_char == 0) {
                    variable_char = (char) ch;
                } else if (ch != variable_char) {
                    System.err.print("Error: multivariable equations are currently unsupported");
                    return false;
                }
                if (curr_fract_part == TOP) {
                    if (this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient == 0) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient = 1;
                    }
                    this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 1;
                } else {
                    if (this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient == 0) {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient = 1;
                    }
                    this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent = 1;
                }
                curr_state = MONOMIAL;
                curr_num = COEFFICIENT;
            } else if (ch == '^') {
                if (!((input_str.charAt(i - 1) >= 'a' && input_str.charAt(i - 1) <= 'z') || (input_str.charAt(i - 1) >= 'A' && input_str.charAt(i - 1) <= 'Z'))) {
                    System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i + 1);
                    return false;
                }
                if (curr_state == MONOMIAL || curr_state == POLYNOMIAL) {
                    if (curr_num != COEFFICIENT) {
                        System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i);
                        return false;
                    }
                    if (curr_fract_part == TOP) {
                        this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                    } else {
                        this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                    }
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }
                    curr_num = EXPONENT;
                }
            } else if (ch == '.') {
                System.out.print("Error: currently does not support decimals\n");
                return false;
            } else if (ch != ' ') {
                for (int j = 0; j < i; j++) {
                    System.out.print(" ");
                }
                System.out.print("^\n");
                System.out.print("Error: invalid or unexpected character at position");
                return false;
            }
        }
        if (is_negative) {
            if (curr_fract_part == TOP) {
                if (curr_num == COEFFICIENT) {
                    this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                } else {
                    this.left_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
                }
            } else if (curr_num == COEFFICIENT) {
                this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
            } else {
                this.left_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
            }
            is_negative = false;
        }
        if (has_equals) {
            curr_monomial = 0;
            curr_polynomial = 0;
            curr_poly_fract = 0;
            curr_state = EXPRESSION;
            curr_num = COEFFICIENT;
            curr_fract_part = TOP;
            for (int i = input_str.indexOf('=') + 1; i < length; i++) {
                ch = input_str.charAt(i);
                if (ch == '(') {
                    curr_monomial = 0;
                    if (curr_state == EXPRESSION || curr_state == POLY_FRACT) {
                        for (int j = i + 1; j < length; j++) {
                            if (input_str.charAt(j) == ')') {
                                break;
                            }
                        }
                    } else {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.add_element();
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.add_element();
                        }
                        curr_polynomial++;
                    }
                    curr_state = POLYNOMIAL;
                } else if (ch == ')') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }
                    if (curr_state == POLY_N_D) {
                        continue;
                    }
                    curr_state = POLY_N_D;
                    curr_monomial = 0;
                } else if (ch == '+') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }
                    if (curr_state == MONOMIAL) {
                        curr_num = COEFFICIENT;
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                        }
                        curr_monomial++;
                    } else if (curr_state == POLY_N_D) {
                        this.right_side.add_element();
                        curr_poly_fract++;
                        curr_num = COEFFICIENT;
                        curr_polynomial = 0;
                        curr_fract_part = TOP;
                        curr_monomial = 0;
                    }
                } else if (ch == '-') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                        is_negative = false;
                    }

                    if (curr_state == MONOMIAL) {
                        curr_num = COEFFICIENT;
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].add_element();
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].add_element();
                        }
                        curr_monomial++;
                        is_negative = true;
                    } else if (curr_state == POLYNOMIAL) {
                        is_negative = true;
                    } else if (curr_state == POLY_N_D) {
                        this.right_side.add_element();
                        curr_poly_fract++;
                        curr_num = COEFFICIENT;
                        curr_polynomial = 0;
                        curr_fract_part = TOP;
                        curr_monomial = 0;
                        is_negative = true;
                    }
                } else if (ch == '/') {
                    if (is_negative) {
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                        }
                    }
                    if (curr_state == MONOMIAL || curr_state == POLY_N_D) {
                        curr_state = POLY_FRACT;
                    }
                    curr_monomial = 0;
                    curr_polynomial = 0;
                    curr_fract_part = BOTTOM;
                    this.right_side.elements[curr_poly_fract].denominator.elements[0].elements[0].coefficient = 0;
                    is_negative = false;
                } else if (ch >= '0' && ch <= '9') {
                    curr_state = MONOMIAL;
                    if (curr_fract_part == TOP) {
                        if (curr_num == COEFFICIENT) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                        } else {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                        }
                    } else {
                        if (curr_num == COEFFICIENT) {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= 10;
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient += ch - '0';
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= 10;
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent += ch - '0';
                        }
                    }
                } else if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')) {
                    if (ch >= 'a' && ch <= 'z') {
                        ch += 'A' - 'a';
                    }
                    if (variable_char == 0) {
                        variable_char = (char) ch;
                    } else if (ch != variable_char) {
                        System.err.print("Error: multivariable equations are currently unsupported");
                        return false;
                    }
                    if (this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient == 0) {
                        this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient = 1;
                    }
                    this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 1;
                    curr_state = MONOMIAL;
                    curr_num = COEFFICIENT;
                } else if (ch == '^') {
                    if (!((input_str.charAt(i - 1) >= 'a' && input_str.charAt(i - 1) <= 'z') || (input_str.charAt(i - 1) >= 'A' && input_str.charAt(i - 1) <= 'Z'))) {
                        System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i + 1);
                        return false;
                    }
                    if (curr_state == MONOMIAL || curr_state == POLYNOMIAL) {
                        if (curr_num != COEFFICIENT) {
                            System.err.printf("Error: invalid or unexpected character '%c' at position %d", ch, i);
                            return false;
                        }
                        if (curr_fract_part == TOP) {
                            this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                        } else {
                            this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent = 0;
                        }
                        if (is_negative) {
                            if (curr_fract_part == TOP) {
                                this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                            } else {
                                this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                            }
                            is_negative = false;
                        }
                        curr_num = EXPONENT;
                    }
                } else if (ch == ' ') {

                } else if (ch == '.') {
                    System.out.print("Error: currently does not support decimals\n");
                    return false;
                } else {
                    for (int j = 0; j < i; j++) {
                        System.out.print(" ");
                    }
                    System.out.print("^\n");
                    System.out.print("Error: invalid or unexpected character at position");
                    return false;
                }
            }
            if (is_negative) {
                if (curr_fract_part == TOP) {
                    if (curr_num == COEFFICIENT) {
                        this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                    } else {
                        this.right_side.elements[curr_poly_fract].numerator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
                    }
                } else if (curr_num == COEFFICIENT) {
                    this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].coefficient *= -1;
                } else {
                    this.right_side.elements[curr_poly_fract].denominator.elements[curr_polynomial].elements[curr_monomial].exponent *= -1;
                }
            }
        }
        System.out.println("done parsing, got: " + this.toString());
        return true;
    }

    //public boolean parse_input() {
    //    return parse_input(userInput.getText());
    //}
    public boolean parse_debug_input(String input) {
        return this.parse_input(input);
    }

    public void set_equal_to_zero() {
        if (right_side.elements[0].numerator.elements[0].elements[0].coefficient == 0) {
            return;
        }
        for (int i = 0; i < right_side.num_elements; i++) {
            right_side.elements[i].numerator.elements[0].elements[0].coefficient *= -1;
            left_side.add_element(right_side.elements[i]);
        }
    }

    public void add_elements() {
        while (left_side.num_elements > 1) {
            left_side.elements[0] = left_side.elements[0].add(left_side.elements[1]);
            left_side.remove_element(1);
        }
    }

    @Override
    public String toString() {
        String str = left_side.toString();
        if (right_side.elements[0].numerator.elements[0].elements[0].coefficient != 0) {
            str = str.concat(" = ");
            str = str.concat(right_side.toString());
        }
        return str;
    }
}

public class EquationSolverUI extends javax.swing.JFrame {

    /**
     * Creates new form EquationSolverUI
     */
    public EquationSolverUI() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        formattingHelp = new javax.swing.JDialog();
        jPanel1 = new javax.swing.JPanel();
        userInput = new javax.swing.JTextField();
        inputLabel = new javax.swing.JLabel();
        btnSolve = new javax.swing.JButton();
        outputFactored = new javax.swing.JTextField();
        outputY_int = new javax.swing.JTextField();
        outputZeros = new javax.swing.JTextField();
        outputHoles = new javax.swing.JTextField();
        outputHasymptotes = new javax.swing.JTextField();
        outputVasymptotes = new javax.swing.JTextField();
        filler1 = new javax.swing.Box.Filler(new java.awt.Dimension(0, 0), new java.awt.Dimension(0, 0), new java.awt.Dimension(0, 32767));
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        jMenuItem1 = new javax.swing.JMenuItem();
        jMenu3 = new javax.swing.JMenu();
        menuFormatting = new javax.swing.JMenuItem();
        menuAbout = new javax.swing.JMenuItem();

        javax.swing.GroupLayout formattingHelpLayout = new javax.swing.GroupLayout(formattingHelp.getContentPane());
        formattingHelp.getContentPane().setLayout(formattingHelpLayout);
        formattingHelpLayout.setHorizontalGroup(
            formattingHelpLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        formattingHelpLayout.setVerticalGroup(
            formattingHelpLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        inputLabel.setLabelFor(inputLabel);
        inputLabel.setText("Enter an equation:");

        btnSolve.setText("Solve");
        btnSolve.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSolveActionPerformed(evt);
            }
        });

        outputFactored.setEditable(false);
        outputFactored.setText("                 ");

        outputY_int.setEditable(false);
        outputY_int.setText("                 ");

        outputZeros.setEditable(false);
        outputZeros.setText("             ");

        outputHoles.setEditable(false);

        outputHasymptotes.setEditable(false);

        outputVasymptotes.setEditable(false);

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(254, 254, 254)
                        .addComponent(btnSolve, javax.swing.GroupLayout.PREFERRED_SIZE, 123, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(outputHoles, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(jPanel1Layout.createSequentialGroup()
                                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(inputLabel)
                                    .addComponent(filler1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(userInput, javax.swing.GroupLayout.DEFAULT_SIZE, 544, Short.MAX_VALUE))
                            .addComponent(outputFactored, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(outputY_int)
                            .addComponent(outputZeros)))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(outputHasymptotes))
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(outputVasymptotes)))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGap(22, 22, 22)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(userInput, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(inputLabel))
                    .addComponent(filler1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(outputFactored, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputY_int, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputZeros, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputHoles, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputHasymptotes, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputVasymptotes, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 101, Short.MAX_VALUE)
                .addComponent(btnSolve)
                .addGap(34, 34, 34))
        );

        jMenu1.setText("File");

        jMenuItem1.setText("Exit");
        jMenuItem1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem1ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem1);

        jMenuBar1.add(jMenu1);

        jMenu3.setText("Help");

        menuFormatting.setText("Formatting");
        menuFormatting.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                menuFormattingActionPerformed(evt);
            }
        });
        jMenu3.add(menuFormatting);

        menuAbout.setText("About");
        jMenu3.add(menuAbout);

        jMenuBar1.add(jMenu3);

        setJMenuBar(jMenuBar1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jMenuItem1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem1ActionPerformed
        System.exit(0);
    }//GEN-LAST:event_jMenuItem1ActionPerformed

    private void btnSolveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSolveActionPerformed
        equation input = new equation();
        input.parse_input(userInput.getText());
        input.set_equal_to_zero();
        input.add_elements();
        input.left_side.factor();
        outputFactored.setText(input.toString());
        Complex[][] solutions = input.left_side.elements[0].find_solutions();
        if (input.left_side.elements[0].denominator.val(0.0) != 0) {
            outputY_int.setText("y-intercept = ".concat(((Double) input.left_side.val(0)).toString()));
        } else {
            outputY_int.setText("no defined y-intercept");
        }
        if (solutions[0].length > 0) {
            String str = "zeros: x = ";
            for (Complex solution : solutions[0]) {
                str = str.concat(solution.toString().concat(", "));
            }
            outputZeros.setText(str);
        }
        if (solutions[1].length > 0) {
            String str = "holes: x=";
            for (Complex solution : solutions[1]) {
                str = str.concat(solution.toString().concat(", "));
            }
            outputHoles.setText(str);
        } else {
            outputHoles.setText("no removable discontinuities (holes)");
        }
        if (solutions[2].length > 0) {
            outputHasymptotes.setText("horizontal asymptote at y = ".concat(solutions[2][0].toString()));
        } else {
            outputHasymptotes.setText("no horizontal asymptote");
        }
        if (solutions[3].length > 0) {
            String str = "vertical asymptote(s) at: x = ";
            for (Complex solution : solutions[3]) {
                str = str.concat(solution.toString().concat(", "));
            }
            outputVasymptotes.setText(str);
        } else {
            outputVasymptotes.setText("no infinite discontinuities (vertical asymptotes)");
        }
    }//GEN-LAST:event_btnSolveActionPerformed

    private void menuFormattingActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuFormattingActionPerformed
    }//GEN-LAST:event_menuFormattingActionPerformed

    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(EquationSolverUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                new EquationSolverUI().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnSolve;
    private javax.swing.Box.Filler filler1;
    private javax.swing.JDialog formattingHelp;
    private javax.swing.JLabel inputLabel;
    private javax.swing.JMenu jMenu1;
    private javax.swing.JMenu jMenu3;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JMenuItem jMenuItem1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JMenuItem menuAbout;
    private javax.swing.JMenuItem menuFormatting;
    private javax.swing.JTextField outputFactored;
    private javax.swing.JTextField outputHasymptotes;
    private javax.swing.JTextField outputHoles;
    private javax.swing.JTextField outputVasymptotes;
    private javax.swing.JTextField outputY_int;
    private javax.swing.JTextField outputZeros;
    private javax.swing.JTextField userInput;
    // End of variables declaration//GEN-END:variables
}
