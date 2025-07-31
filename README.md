# ECDLP_Pollard_Rho

# Elementary Pollard's Rho for ECDLP
This repository contains a basic, elementary implementation of Pollard's Rho algorithm for solving the **Elliptic Curve Discrete Logarithm Problem (ECDLP)**. The code is designed to be simple and clear for educational purposes.
-----
## ⚠️ **Security Warning: For Educational Use Only** ⚠️
This code is intended for **learning and demonstration purposes only**. It is a simplified, "textbook" implementation and is **NOT CRYPTOGRAPHICALLY SECURE**. Real-world elliptic curve cryptography relies on the difficulty of solving the ECDLP; this code is an attack on that problem and demonstrates how it can be broken for small parameters.
**DO NOT use this code in any production environment, for real-world applications, or with sensitive data.**
Key limitations include:
  * It is not optimized and is only feasible for curves with a **very small** group order.
  * It lacks necessary safeguards and may fail or crash on certain inputs.
  * It is not a robust or reliable cryptographic tool.
-----
## Purpose
The primary goal of this repository is to provide a straightforward example of how Pollard's Rho algorithm can be applied to solve the ECDLP. It is a learning tool for students and enthusiasts interested in elliptic curve cryptography and number theory algorithms.
-----
## Usage
The elliptic curve parameters (such as the prime field, curve coefficients, and base point order) and the points $P$ and $Q$ must be defined within the script itself.
Once the parameters are configured in the code, you can execute the script.

