from numpy import max, argmax


class problem:
    def __init__(self, L, W1, W2, x):
        """
        L: The length of the beam in metres.
        W1: The first load in kN
        W2: The second load in kN
        x: The spacing between the loads in metres.
        """

        self.L = L
        self.W1, self.W2 = W1 * 1000, W2 * 1000  # converted kN to N
        self.x = x
        self.valid_problem = (
            x < L
        )  # Spacing between the loads shouldn't exceed the beam length
        self.x_by_L = x / L
        self.U1 = self.W1 / (self.W1 + self.W2)
        self.U2 = self.W2 / (self.W1 + self.W2)

    def BM_infline(self, p, a):
        """
        p: Distance from support A where BM is being evaluated
        a: Distance of the load from support A

        Returns:
        The value of the influence line for the bending moment.
        """
        if a < 0 or a > self.L:
            return 0
        elif a <= p:
            return a * (1 - p / self.L)
        else:
            return (1 - a / self.L) * p

    def SF_infline(self, p, a):
        """
        p: Distance from support A where BM is being evaluated
        a: Distance of the load from support A

        Returns:
        The value of the influence line for the Shear Force.
        """
        if a < 0 or a > self.L:
            return 0
        elif a < p:
            return -(a / self.L)
        else:
            return 1 - a / self.L

    def SF_calc(self, p, a):
        a1, a2 = a, a + self.x  # distance of both loads from A

        # Total SF = sum (load * unit contribution of load)
        return self.W1 * self.SF_infline(p, a1) + self.W2 * self.SF_infline(p, a2)

    def BM_calc(self, p, a):
        a1, a2 = a, a + self.x

        # Total BM = sum (load * unit contribution of load)
        return self.W1 * self.BM_infline(p, a1) + self.W2 * self.BM_infline(p, a2)

    def get_max_A(self):  # calculate maximum reaction at A
        return max([self.W2, self.W1 + self.W2 * (1 - self.x_by_L)])

    def get_max_B(self):  # calculate maximum reaction at B
        return max([self.W1, self.W2 + self.W1 * (1 - self.x_by_L)])

    def get_BM_01(self):
        return self.W2 * self.x * (1 - self.x_by_L)

    def get_SF_01(self):
        return 0.5 * (self.W1 + self.W2 + self.x_by_L * abs(self.W1 - self.W2))

    def get_max_SF(self):
        if self.W1 > self.W2:
            return 0, self.W1 + self.W2 * (1 - self.x_by_L)
        else:
            return self.L, self.W2 + self.W1 * (1 - self.x_by_L)

    def get_max_BM(self):
        # Critical points enumerated in 6.6
        crit_p = [
            [
                0.5 * (self.L + self.U1 * self.x),
                0.5 * (self.L + self.U1 * self.x) - self.x,
            ],
            [0.5 * (self.L - self.U2 * self.x), 0.5 * (self.L - self.U2 * self.x)],
            [self.L - self.x, self.L - self.x],
            [self.x, 0],
        ]

        # Calculate the BM values by passing to BM_calc
        BM_values = [self.BM_calc(x[0], x[1]) for x in crit_p]

        # Maximize
        max_BM = max(BM_values)
        max_BM_pos = crit_p[argmax(BM_values)][0]

        # Return the max value and where it occurs.
        return max_BM_pos, max_BM


if __name__ == "__main__":
    # L = input("Enter the length of the Beam (m): ")
    # W1 = input("Enter the first load (kN): ")
    # W2 = input("Enter the second load (kN): ")
    # x = input("Enter the spacing x (m): ")

    # Toy problem
    osdag_task = problem(10, 3, 1, 4)

    # Maximum Reaction at A
    max_A = osdag_task.get_max_A()
    print(f"Maximum reaction at support A: {max_A} N")

    # Maximum Reaction at B
    max_B = osdag_task.get_max_B()
    print(f"Maximum reaction at support B: {max_B} N")

    # Get BM_01
    BM_01 = osdag_task.get_BM_01()
    print(f"Value of BM_01: {BM_01} N.m")

    # Get SF_01
    SF_01 = osdag_task.get_SF_01()
    print(f"Value of SF_01: {SF_01} N")

    # Maximum value of Shear Force
    pos_SF, max_SF = osdag_task.get_max_SF()
    print(f"Max Shear force {max_SF} N occurs at {pos_SF} m from support A.")

    # Maximum value of Bending Moment
    pos_BM, max_BM = osdag_task.get_max_BM()
    print(f"Max Bending moment {max_BM} occurs at {pos_BM} m from support A.")
