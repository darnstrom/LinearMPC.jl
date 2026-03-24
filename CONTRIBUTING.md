# Contributing to LinearMPC.jl

First off, thank you for taking the time to contribute! It is people like you who make the Julia control systems ecosystem great.

All types of contributions are welcome: from reporting bugs and improving documentation to submitting performance enhancements and new MPC formulations.

---

## How Can I Contribute?

### Reporting Bugs
Before creating a bug report, please check the [Issues](https://github.com/darnstrom/LinearMPC.jl/issues) tab to see if the problem has already been reported.

When filing an issue, please include:
* **Julia version** and **Package version**.
* A **Minimal Working Example (MWE)**: A small, self-contained script that reproduces the error.
* The full error trace (if applicable).

### Suggesting Enhancements
If you have an idea for a new feature (e.g., a new solver interface or support for time-varying systems):
1.  Open an issue to discuss it first.
2.  Provide a clear description of the use case and how it benefits the package.

### Pull Requests (PRs)
1.  **Fork** the repository and create your branch from `main`.
2.  **Install dependencies**: Run `import Pkg; Pkg.activate("."); Pkg.instantiate()`.
3.  **Implement changes**: Ensure your code follows the existing style (concise and type-stable).
4.  **Add Tests**: If you add a feature, add a corresponding test file in `/test`.
5.  **Run Tests**: Ensure everything passes by running:
    ```julia
    import Pkg; Pkg.test("LinearMPC")
    ```
6.  **Submit**: Open a PR with a concise title and a description of your changes.

---

## Coding Standards

To keep the codebase maintainable, please keep the following in mind:

* **Type Stability**: Since MPC is often time-critical, ensure that inner loops are type-stable and avoid unnecessary allocations.
* **Docstrings**: Use Julia's native docstring format for all exported functions. Include arguments, return types, and a small example.

---

## Code of Conduct
By participating in this project, you agree to maintain a respectful, inclusive, and professional environment. 
