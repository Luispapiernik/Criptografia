ABC = [' ', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
       'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']


def modulo(function):
    def wrapper(self, number):
        result = function(self, number)
        return Zn(result, self.base)
    return wrapper


def decorateClass(cls):
    functions = ['__add__', '__sub__', '__mul__', '__pow__',
                 '__radd__', '__rsub__', '__rmul__']

    for function in functions:
        setattr(cls, function, modulo(getattr(cls, function)))

    return cls


def gcd(a, b):
    if a % b == 0:
        return b
    return gcd(b, a % b)


class MathError(Exception):
    pass


class NotInvertibleError(MathError):
    def __init__(self):
        self.message = 'Element no invertible'
        super(NotInvertibleError, self).__init__(self.message)


class DimensionError(MathError):
    def __init__(self):
        self.message = 'The matrices have incompatible dimension'
        super(DimensionError, self).__init__(self.message)


@decorateClass
class Zn(int):
    def __new__(cls, value, base=26):
        return super(Zn, cls).__new__(cls, value % base)

    def __init__(self, value, base=26):
        self.base = base
        self.isinvertible = 1 == gcd(self, self.base)

    def __repr__(self):
        return '(%d, %d)' % (self, self.base)

    def inverse(self):
        if not self.isinvertible:
            raise NotInvertibleError

        inv_1 = 0
        inv_2 = 1

        number = self
        base = self.base

        while base % number != 0:
            inv_1, inv_2 = inv_2, inv_1 - inv_2 * (base // number)
            number, base = base % number, number

        return Zn(inv_2, self.base)

    def __truediv__(self, divisor):
        if isinstance(divisor, int):
            divisor = Zn(divisor, self.base)
        return self * divisor.inverse()

    def __rtruediv__(self, dividendo):
        return dividendo * self.inverse()

    def __pow__(self, exponent):
        previus_power = int(self)
        if exponent < 0:
            exponent = abs(exponent)
            previus_power = int(self.inverse())

        c = 1
        new_exponent = 0
        while exponent:
            if exponent & 1:
                previus_power = (
                    previus_power ** (2 ** new_exponent)) % self.base
                c *= previus_power
                c %= self.base
                new_exponent = 0

            new_exponent += 1

            exponent >>= 1

        return c


class MZn(object):
    def __init__(self, matrix, base=26, rows=0, columns=0, modulo=True):
        if modulo:
            matrix = [[Zn(number, base) for number in row] for row in matrix]

        self.matrix = matrix
        self.rows = rows or len(matrix)
        self.columns = columns or len(matrix[0])

        self.base = base
        digits = len(str(base))
        self.parser = ('{:%dd}' % (digits + 1)) * self.columns

    def __repr__(self):
        return str(self.matrix)

    def __str__(self):
        string = ''
        for line in self.matrix:
            string += self.parser.format(*line)
            string += '\n'

        return string[:-1]

    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            row, column = key
            if isinstance(row, int) and isinstance(column, int):
                return self.matrix[row][column]
        return MZn(self.matrix[key], base=self.base, modulo=False,
                   columns=self.columns)

    def transpose(self):
        transposed = [[0] * self.rows for i in range(self.columns)]

        for i in range(self.rows):
            for j in range(self.columns):
                transposed[j][i] = self.matrix[i][j]

        return MZn(transposed, base=self.base, rows=self.columns,
                   columns=self.rows, modulo=False)

    def zip(self, matrix, binaryFunction):
        if matrix.rows < self.rows or matrix.columns < self.columns:
            raise DimensionError

        result = [[0] * self.columns for i in range(self.rows)]

        for i in range(self.rows):
            for j in range(self.columns):
                result[i][j] = binaryFunction(self[i, j], matrix[i, j])

        return MZn(result, base=self.base, rows=self.rows,
                   columns=self.columns, modulo=False)

    def __add__(self, matrix):
        if matrix.rows == self.rows and matrix.columns == self.columns:
            return self.zip(matrix, lambda x, y: x + y)

        raise DimensionError

    def __sub__(self, matrix):
        if matrix.rows == self.rows and matrix.columns == self.columns:
            return self.zip(matrix, lambda x, y: x - y)

        raise DimensionError

    def __mul__(self, matrix):
        if not isinstance(matrix, MZn):
            matrix = MZn(matrix, base=self.base)

        if self.columns != matrix.rows:
            raise DimensionError

        product = [[0] * matrix.columns for i in range(self.rows)]

        for i in range(self.rows):
            for j in range(matrix.columns):
                component_ij = 0
                for k in range(self.columns):
                    component_ij += self.matrix[i][k] * matrix[k, j]

                product[i][j] = component_ij

        return MZn(product, base=self.base, rows=self.rows,
                   columns=matrix.columns, modulo=False)

    def __rmul__(self, matrix):
        matrix = MZn(matrix, base=self.base)

        return matrix * self

    def matMul(self, matrix, transpose=False):
        if not isinstance(matrix, MZn):
            matrix = MZn(matrix, base=self.base)

        if transpose:
            matrix = matrix.transpose()

        return self * matrix

    def dot(self, matrix):
        if matrix.rows == self.rows and matrix.columns == self.columns:
            return self.zip(matrix, lambda x, y: x * y)

        raise DimensionError

    def inverse(self):
        raise NotImplementedError

    def identity(self):
        one = Zn(1, self.base)
        zero = Zn(0, self.base)
        result = [[zero] * i + [one] + [zero] * (self.columns - i - 1)
                  for i in range(self.rows)]

        return MZn(result, base=self.base, rows=self.rows,
                   columns=self.columns)

    def __pow__(self, exponent):
        raise NotImplementedError
        if exponent < 0:
            pass

        result = self.identity()

        for i in range(exponent):
            result *= self

        return result


def testInverse():
    a = Zn(3)
    b = Zn(8)
    c = Zn(11)

    ai = a ** -1
    bi = b ** -1
    ci = c ** -1

    print(a, b, c)
    print(ai, bi, ci)
    print(a * ai, b * bi, c * ci)

    print(a ** 4, b ** 2, 11 ** 1)


def testOperations():
    a = Zn(7)

    print(a, a * 7, 7 * a)
    print(a, a ** 2, 2 ** a)


def main():
    testInverse()
    testOperations()


if __name__ == '__main__':
    main()
