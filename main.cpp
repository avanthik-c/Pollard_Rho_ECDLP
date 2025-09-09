#include <iostream>
#include <stdexcept> // Required for std::runtime_error
#include <cmath>     // Required for std::abs

// Note: For real-world cryptography, a BigInt (Arbitrary Precision Arithmetic)
// library is essential because 'long long' is too small. This implementation uses 'long long'
// for educational purposes to demonstrate the algorithm's structure.
using int_type = long long;

// Utility function for modular arithmetic to handle potential negative results.
int_type mod(int_type a, int_type n) {
    if (n == 0) throw std::runtime_error("Modulo by zero");
    return (a % n + n) % n;
}

// Extended Euclidean Algorithm to find modular inverse for integers.
// Returns x such that (a*x) % m = 1.
// Throws an error if no inverse exists.
int_type integerModInverse(int_type a, int_type m) {
    int_type m0 = m, t, q;
    int_type x0 = 0, x1 = 1;

    if (m == 1) return 0;
    
    a = mod(a, m);

    while (a > 1) {
        if (m == 0) throw std::runtime_error("No modular inverse exists (division by zero).");
        q = a / m;
        t = m;
        m = a % m, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }

    if (x1 < 0) x1 += m0;

    return x1;
}

/**
 * @class Point
 * Represents a point (x, y) on an elliptic curve.
 *
 * @oop_concept Encapsulation
 * This struct encapsulates the data related to a point on a curve (its coordinates x, y,
 * and its state as the point at infinity) into a single unit.
 *
 * @oop_concept Abstraction
 * It abstracts the concept of a point, hiding the details of its representation. Users
 * interact with it as a single entity rather than separate coordinate variables.
 */
struct Point {
    int_type x;
    int_type y;
    bool is_infinity;

    // Default constructor for the point at infinity
    Point() : x(0), y(0), is_infinity(true) {}

    // Constructor for a regular point on the curve
    Point(int_type x_val, int_type y_val) : x(x_val), y(y_val), is_infinity(false) {}

    /**
     * @oop_concept Polymorphism (Ad-hoc / Overloading)
     * We overload the '==' operator to provide a specific implementation for comparing
     * two Point objects. This allows for intuitive comparison syntax (e.g., `if (P == Q)`).
     */
    bool operator==(const Point& other) const {
        if (is_infinity && other.is_infinity) {
            return true;
        }
        if (is_infinity || other.is_infinity) {
            return false;
        }
        return (x == other.x && y == other.y);
    }
};

/**
 * @class EllipticCurve
 * Represents an elliptic curve defined by y^2 = x^3 + ax + b (mod p) over a prime field.
 *
 * @oop_concept Encapsulation
 * This class encapsulates the curve's parameters (a, b, p) and the specialized
 * mathematical operations (addition, doubling) that apply to points on this curve.
 *
 * @oop_concept Abstraction
 * It hides the complex formulas of elliptic curve arithmetic behind a simple
 * interface (`add`, `doublePoint`). The user does not need to know the underlying
 * mathematics to use the class.
 */
class EllipticCurve {
private:
    int_type a;
    int_type p; // Prime modulus of the field

public:
    EllipticCurve(int_type param_a, int_type prime_mod) : a(param_a), p(prime_mod) {}

    // Doubles a point P on the curve (calculates 2P).
    Point doublePoint(Point pt) {
        if (pt.is_infinity || pt.y == 0) {
            return Point(); // Point at infinity
        }

        // Calculate slope of the tangent line (m)
        // m = (3*x^2 + a) * (2*y)^-1 mod p
        int_type numerator = mod(3 * pt.x * pt.x + a, p);
        int_type denominator = mod(2 * pt.y, p);
        int_type inv_denominator = integerModInverse(denominator, p);

        int_type m = mod(numerator * inv_denominator, p);

        // Calculate coordinates of the new point (x_r, y_r)
        // x_r = m^2 - 2*x mod p
        // y_r = m*(x - x_r) - y mod p
        int_type xr = mod(m * m - 2 * pt.x, p);
        int_type yr = mod(m * (pt.x - xr) - pt.y, p);

        return Point(xr, yr);
    }

    // Adds two distinct points P and Q on the curve.
    Point add(Point pt1, Point pt2) {
        if (pt1.is_infinity) return pt2;
        if (pt2.is_infinity) return pt1;

        if (pt1 == pt2) {
            return doublePoint(pt1);
        }
        
        // If P.x == Q.x but P.y != Q.y (i.e., P.y = -Q.y mod p), the line is vertical.
        if (pt1.x == pt2.x) {
            return Point(); // Result is the point at infinity
        }

        // Calculate slope of the line between P and Q (m)
        // m = (pt2.y - pt1.y) * (pt2.x - pt1.x)^-1 mod p
        int_type numerator = mod(pt2.y - pt1.y, p);
        int_type denominator = mod(pt2.x - pt1.x, p);
        int_type inv_denominator = integerModInverse(denominator, p);
        
        int_type m = mod(numerator * inv_denominator, p);

        // Calculate coordinates of the new point (x_r, y_r)
        // x_r = m^2 - pt1.x - pt2.x mod p
        // y_r = m*(pt1.x - x_r) - pt1.y mod p
        int_type xr = mod(m * m - pt1.x - pt2.x, p);
        int_type yr = mod(m * (pt1.x - xr) - pt1.y, p);

        return Point(xr, yr);
    }
};

/**
 * @class ECDLP_Solver
 * Abstract base class defining the interface for any ECDLP solving algorithm.
 *
 * @oop_concept Inheritance (Base Class)
 * This class serves as a blueprint for concrete ECDLP solver classes.
 *
 * @oop_concept Abstraction (Abstract Class)
 * It cannot be instantiated. It defines a contract (`solve` method) that
 * any derived class must implement.
 */
class ECDLP_Solver {
public:
    virtual ~ECDLP_Solver() = default;
    
    /**
     * @oop_concept Polymorphism (Runtime via Pure Virtual Function)
     * This pure virtual function allows different ECDLP algorithms to be used
     * interchangeably through a base class pointer. The correct `solve` method
     * is called at runtime based on the actual object's type.
     */
    virtual int_type solve() = 0;
};

/**
 * @class PollardRhoECDLPSolver
 * Implements the Pollard's Rho algorithm to solve the ECDLP.
 *
 * @oop_concept Inheritance (Derived Class)
 * This class inherits from the abstract `ECDLP_Solver` class, providing a concrete
 * implementation of the ECDLP solving algorithm.
 *
 * @oop_concept Composition ("has-a" relationship)
 * This class contains an `EllipticCurve` object. It uses this
 * object to perform its work, delegating curve arithmetic to it.
 */
class PollardRhoECDLPSolver : public ECDLP_Solver {
private:
    EllipticCurve curve;
    Point P;       // Base point
    Point Q;       // Target point (Q = kP)
    int_type order; // Order of the subgroup generated by P

public:
    PollardRhoECDLPSolver(int_type a, int_type p, Point base_P, Point target_Q, int_type subgroup_order)
        : curve(a, p), P(base_P), Q(target_Q), order(subgroup_order) {}

    // Defines the pseudo-random walk for Pollard's Rho.
    // X_{i+1} is computed based on X_i.
    // We also update coefficients a, b such that X_i = aP + bQ.
    void advance(Point& X, int_type& a, int_type& b) {
        // Partition points into 3 sets based on x-coordinate for pseudo-randomness
        switch (X.x % 3) {
            case 0: // Add P
                X = curve.add(X, P);
                a = mod(a + 1, order);
                break;
            case 1: // Double
                X = curve.doublePoint(X);
                a = mod(a * 2, order);
                b = mod(b * 2, order);
                break;
            case 2: // Add Q
                X = curve.add(X, Q);
                b = mod(b + 1, order);
                break;
        }
    }

    /**
     * @oop_concept Polymorphism (Runtime via Overriding)
     * This method overrides the pure virtual function from the base `ECDLP_Solver` class.
     */
    int_type solve() override {
        // Initial state for tortoise and hare
        Point tortoise_X = P;
        int_type tortoise_a = 1, tortoise_b = 0;

        Point hare_X = P;
        int_type hare_a = 1, hare_b = 0;

        // First step for hare
        advance(hare_X, hare_a, hare_b);

        int max_iter = 2 * order; // Should find a cycle within the order
        for (int i = 0; i < max_iter; ++i) {
            if (tortoise_X == hare_X) {
                // Collision found! X_i = X_{2i}
                // a_i*P + b_i*Q = a_{2i}*P + b_{2i}*Q
                // (b_i - b_{2i})*Q = (a_{2i} - a_i)*P
                // Since Q = kP, we have:
                // k*(b_i - b_{2i})*P = (a_{2i} - a_i)*P
                // k*(b_i - b_{2i}) = (a_{2i} - a_i) mod order
                // k = (a_{2i} - a_i) * (b_i - b_{2i})^-1 mod order
                int_type numerator = mod(hare_a - tortoise_a, order);
                int_type denominator = mod(tortoise_b - hare_b, order);
                
                if (denominator == 0) {
                    // This is a failure case, very rare.
                    // It means tortoise_b == hare_b mod order.
                    // Try again with a different random walk.
                    return -1;
                }

                int_type inv_denominator = integerModInverse(denominator, order);
                return mod(numerator * inv_denominator, order);
            }

            // Tortoise moves one step
            advance(tortoise_X, tortoise_a, tortoise_b);
            // Hare moves two steps
            advance(hare_X, hare_a, hare_b);
            advance(hare_X, hare_a, hare_b);
        }
        
        return -1; // No solution found within max iterations
    }
};

int main() {
    // ---- Problem Setup for ECDLP ----
    // Solve for k in Q = kP on the curve y^2 = x^3 + ax + b (mod p)

    // Example: Curve y^2 = x^3 + 2x + 2 (mod 17)
    // The value of 'b' is implicitly defined by the curve parameters.
    int_type p = 17;
    int_type a = 2;
    // Base point P
    Point P(5, 1);
    // Order of the subgroup generated by P is 19.
    int_type order = 19;
    // Target point Q. We want to find k such that Q = kP.
    // For this example, let's use Q = 10P = (10, 6)
    Point Q(10, 6);

    std::cout << "--- Solving ECDLP ---\n";
    std::cout << "Curve: y^2 = x^3 + " << a << "x + b (mod " << p << ")\n";
    std::cout << "Base Point P = (" << P.x << ", " << P.y << ")\n";
    std::cout << "Target Point Q = (" << Q.x << ", " << Q.y << ")\n";
    std::cout << "Subgroup Order = " << order << "\n\n";

    PollardRhoECDLPSolver solver(a, p, P, Q, order);
    
    // Using the abstract base class pointer to demonstrate runtime polymorphism.
    ECDLP_Solver* p_solver = &solver;
    int_type k = p_solver->solve();

    if (k != -1) {
        std::cout << "SUCCESS: Found discrete logarithm k = " << k << std::endl;
        std::cout << "Verification: " << k << "P = Q\n";
    } else {
        std::cout << "FAILURE: Could not find the discrete logarithm." << std::endl;
    }

    return 0;
}

