x = 7
y = 2


class InvalidInputException(Exception):
    pass

try:
    print(f"{x/y}")
except ZeroDivisionError as ex:
    print("No Zeros")
    raise InvalidInputException('Invalid Input')
finally:
    print("lalallalal")