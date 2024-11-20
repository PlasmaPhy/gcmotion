# Higher-order central difference for second derivatives
def higher_order_second_derivative(f, x: float, y: float, dx: float, dy: float, respect_to: str):
    if respect_to == "x":
        return (
            -f(x + 2 * dx, y)
            + 16 * f(x + dx, y)
            - 30 * f(x, y)
            + 16 * f(x - dx, y)
            - f(x - 2 * dx, y)
        ) / (12 * dx**2)
    elif respect_to == "y":
        return (
            -f(x, y + 2 * dy)
            + 16 * f(x, y + dy)
            - 30 * f(x, y)
            + 16 * f(x, y - dy)
            - f(x, y - 2 * dy)
        ) / (12 * dy**2)
    else:  # Mixed derivative
        return (f(x + dx, y + dy) - f(x + dx, y - dy) - f(x - dx, y + dy) + f(x - dx, y - dy)) / (
            4 * dx * dy
        )
