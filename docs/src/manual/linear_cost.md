# Linear Cost

LinearMPC supports a linear cost term $l_k^T u_k$ that penalizes control inputs, enabling economic MPC where control actions have time-varying costs. A typical application is heating or cooling systems where electricity prices vary throughout the day.

## Example: Heat Pump with Time-Varying Electricity Prices

This example demonstrates using linear cost to optimize heating while accounting for electricity price variations.

### Price Generation

First, we define a function to generate realistic electricity prices with daily patterns:

```@tabsetup linear_cost
# julia
function generate_test_prices(n_days=7)
    n_per_day = 24
    n = n_days * n_per_day
    prices = zeros(n)
    for day in 1:n_days
        t = range(0, 24, length=n_per_day)
        daily_factor = 0.8 + 0.4 * rand()
        base = 0.5 .+ 0.2 * sin.(2π * (t .- 6) / 24)
        morning_shift = 8 + randn() * 0.5
        evening_shift = 18 + randn() * 0.5
        morning_peak = (0.4 + 0.2 * rand()) * exp.(-((t .- morning_shift) .^ 2) / 2)
        evening_peak = (0.3 + 0.2 * rand()) * exp.(-((t .- evening_shift) .^ 2) / 3)
        night_dip = -0.1 * exp.(-((t .- 3) .^ 2) / 4)
        day_prices = daily_factor * (base + morning_peak + evening_peak + night_dip)
        day_prices = day_prices .+ 0.05 * randn(n_per_day)
        idx_start = (day - 1) * n_per_day + 1
        idx_end = day * n_per_day
        prices[idx_start:idx_end] = day_prices
    end
    prices = max.(prices, 0.1)
    prices = 0.2 .+ 1.5 * (prices .- minimum(prices)) / (maximum(prices) - minimum(prices))
    return prices
end
# python
import numpy as np

def generate_test_prices(n_days=7):
    n_per_day = 24
    n = n_days * n_per_day
    prices = np.zeros(n)
    for day in range(n_days):
        t = np.linspace(0, 24, n_per_day)
        daily_factor = 0.8 + 0.4 * np.random.rand()
        base = 0.5 + 0.2 * np.sin(2 * np.pi * (t - 6) / 24)
        morning_shift = 8 + np.random.randn() * 0.5
        evening_shift = 18 + np.random.randn() * 0.5
        morning_peak = (0.4 + 0.2 * np.random.rand()) * np.exp(-((t - morning_shift) ** 2) / 2)
        evening_peak = (0.3 + 0.2 * np.random.rand()) * np.exp(-((t - evening_shift) ** 2) / 3)
        night_dip = -0.1 * np.exp(-((t - 3) ** 2) / 4)
        day_prices = daily_factor * (base + morning_peak + evening_peak + night_dip)
        day_prices = day_prices + 0.05 * np.random.randn(n_per_day)
        idx_start = day * n_per_day
        prices[idx_start:idx_start + n_per_day] = day_prices
    prices = np.maximum(prices, 0.1)
    prices = 0.2 + 1.5 * (prices - prices.min()) / (prices.max() - prices.min())
    return prices
```

### Room Temperature Model

We model a simple room heated by a heat pump. The temperature dynamics are first-order:
```math
T_{k+1} = a \cdot T_k + b \cdot u_k
```
where $T$ is room temperature, $u$ is heating power, $a < 1$ represents heat loss, and $b$ is the heating efficiency.

```@tabsetup linear_cost
# julia
using LinearMPC
using Random
Random.seed!(0)

# First-order temperature dynamics (discrete-time, 1 hour sampling)
a = 0.95  # Heat retention
b = 4.0   # Heating efficiency

A = [a;;]
B = [b;;]
C = [1.0;;]

Np = 24  # 24-hour prediction horizon

mpc = LinearMPC.MPC(A, B; C, Np=Np, Nc=Np)
set_bounds!(mpc; umin=[0.0], umax=[1.0])
set_objective!(mpc; Q=[0.05])

# Enable linear cost
mpc.settings.linear_cost = true
setup!(mpc)
# python
from lmpc import MPC
import numpy as np
np.random.seed(0)

# First-order temperature dynamics (discrete-time, 1 hour sampling)
a = 0.95  # Heat retention
b = 4.0   # Heating efficiency

A = np.array([[a]])
B = np.array([[b]])
C = np.array([[1.0]])

Np = 24  # 24-hour prediction horizon

mpc = MPC(A, B, C=C, Np=Np, Nc=Np)
mpc.set_bounds(umin=[0.0], umax=[1.0])
mpc.set_objective(Q=[0.05])

# Enable linear cost
mpc.settings({"linear_cost": True})
mpc.setup()
```

### Simulation

We simulate 3 days of operation, maintaining a temperature setpoint while minimizing electricity cost:

```@tabsetup linear_cost
# julia
N_sim = 72  # 3 days

# Generate electricity prices
prices = generate_test_prices(3)

# Reference temperature
r = fill(21.0, 1, N_sim)

# Linear cost is the electricity price (nu × N_sim)
l = prices'

# Run simulation
sim = Simulation(mpc; x0=[18.0], N=N_sim, r, l)
# python
from lmpc import Simulation

N_sim = 72  # 3 days

# Generate electricity prices
prices = generate_test_prices(3)

# Reference temperature
r = np.full((1, N_sim), 21.0)

# Linear cost is the electricity price (nu × N_sim)
l = prices.reshape(1, -1)

# Run simulation
sim = Simulation(mpc, x0=[18.0], N=N_sim, r=r, l=l)
```

### Results

```@tabexample linear_cost
# julia
using Plots
plot(sim)
# python
import matplotlib.pyplot as plt
plt.plot(sim.ts, sim.ys.T, label="Temperature")
plt.xlabel("Hour"); plt.ylabel("Temperature (°C)")
plt.show()
```

The controller pre-heats when electricity is cheap and reduces heating during expensive peak hours, while keeping the temperature close to the setpoint.

We can also visualize the electricity prices alongside the heating pattern:

```@tabexample linear_cost
# julia
p1 = plot(sim.us', label="Heating power", ylabel="Power")
p2 = plot(l', label="Electricity price", ylabel="Price", xlabel="Hour")
plot(p1, p2, layout=(2,1))
# python
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
ax1.plot(sim.ts, sim.us.T, label="Heating power")
ax1.set_ylabel("Power"); ax1.legend()
ax2.plot(sim.ts, l.T, label="Electricity price")
ax2.set_ylabel("Price"); ax2.set_xlabel("Hour"); ax2.legend()
plt.tight_layout(); plt.show()
```

Notice how heating tends to increase when prices are low (typically at night) and decrease during peak price periods. We can demonstrate this further by plotting the prices vs. the heating power, which should show an inverse relationship:

```@tabexample linear_cost
# julia
scatter(l', sim.us'; xlabel="Electricity Price", ylabel="Heating Power", title="Heating Power vs. Electricity Price")
# python
plt.scatter(l.T, sim.us.T)
plt.xlabel("Electricity Price"); plt.ylabel("Heating Power")
plt.title("Heating Power vs. Electricity Price")
plt.show()
```


### Conclusion
Using a linear cost in the MPC objective allows the controller to optimize heating actions based on time-varying electricity prices, leading to cost-effective operation while maintaining comfort.

The example can be extended in several directions
- Time varying reference temperature (e.g., lower at night)
- Adding constraints on temperature (e.g., minimum temperature)
