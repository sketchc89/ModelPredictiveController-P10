# CarND-Controls-MPC

Self-Driving Car Engineer Nanodegree Program

---

![MPC 70](./Docs/mpc_71.gif)

## Implementation

This project demonstrates Model Predictive Control to control a simulated vehicle to a set trajectory around a simulated track. The model has a state, a set of constraint, and a cost assigned to its state. A non-linear solver is used to satisfy these constraints while reducing the cost.

### State

Following some suggestions on Udacity's forums, I decided to use a simplified state. The state only includes the steering angle and throttle value. The simplified state makes the cost optimization easier to solve. Using the full state has 4 times as many variables and takes longer to compute. The longer time increases the latency of the system which causes the car to control poorly. The other aspects of the state (x, y, psi, cross track error, velocity, psi error) are calculated to compute cost, but they are not included in the state. This state update is run everytime the solver cycles in fg_eval.

The kinematic model uses the following equations:

* x<sub>t+1</sub> = x<sub>t</sub> \* cos(psi<sub>t</sub>) \* dt
* y<sub>t+1</sub> = y<sub>t</sub> \* sin(psi<sub>t</sub>) \* dt
* psi<sub>t+1</sub> = psi<sub>t</sub> + (vel<sub>t</sub> \* dt) / L<sub>f</sub>
* vel<sub>t+1</sub> = vel<sub>t</sub> + acc<sub>t</sub> \* dt
* cte<sub>t+1</sub> = f(x<sub>t</sub>) + y<sub>t</sub> + v<sub>t</sub> \* sin(psi<sub>t</sub>) \* dt
* psi_err<sub>t+1</sub> = psi<sub>t</sub> - psidest<sub>t</sub> + (vel<sub>t</sub> \* dt) / L<sub>f</sub>

where (all metrics are relative to the center of gravity of the car ):

* x = longitudinal position of the car (forward positive)
* y = lateral position of the car (right positive)
* psi = heading position of car (clockwise positive)
* vel = velocity of car (forward positive)
* cte = distance between the lateral position of the trajectory and the lateral position of the car
* psierr = angle between the heading of the proposed trajectory of the car and the heading position of the car
* del = steering angle of car
* acc = acceleration of car
* dt = time between prediction steps
* Lf = distance between center of gravity of the car and the front of the vehicle
* f(x) = the lateral position of the trajectory given the longitudinal position of the car
* psidest = the heading of the proposed trajectory of the car given the cars longitudinal position

### Cost

Cost is added from the following sources:

1. Velocity relative to a reference velocity (70 mph)
2. Cross-track error
3. Heading (psi) error
4. Magnitude of steering angle
5. Magnitude of throttle
6. Change in steering angle
7. Change in throttle
8. Product of magnitude of steering angle and magnitude of throttle

The parameters were weighted relative to one another in order with 0 added latency, then adjusted after adding 100 ms of latency. The most notable cost added to the system is #8. This cost penalizes high velocity sharp turns which causes the car to tend to understeer. Having a small bit of understeer is more desirable than oversteer because oversteering systems are dynamically unstable.

### Timing Parameters

The number of timesteps was originally set to 25 and time between steps was set to 0.05. As the vehicle reference speed was increased, I reduced the number of timesteps to 10 and the time between steps to 0.1. I reduced the number of timesteps in order to reduce computation time of the Ipopt solver. Once I reduced the number of timesteps, I had to increase the time between steps in order to maintain the same look ahead distance. A smaller dt is better because it allows the system to update the vehicles projected path more quickly, but it has to be weighed against the extra computation time that short cycles incur.

### Polynomial fitting

A list of points representing points along the trajectory and the vehicles position and heading were received from the simulation. The rotation matrix is R(psi) = [cos(psi) -sin(psi); sin(psi) cos(psi)]. We would like the trajectory to be fitted to a polynomial relative to the center of gravity of the car (car Cg = 0,0,0; x,y,yaw). In order to do this we rotate the global coordinate system by negative psi because we're rotating from psi back to 0.

Once the global coordinates are converted to vehicle coordinates, we use the trajectory points to solve for a cubic polynomial. The method used is a least square approximation by QR factorization. QR factorization was done with a householder transformation implementation from the Eigen library, because it was used in the cloned scaffolding of this project.

### Latency adjustment

100 ms of latency was added to the system to simulate the effect of sending commands to actuators (this seems excessive). In order to compensate for latency, the first measurement of state was calculated as if it were one step into the future. This takes place in FG_eval::KinematicModel. The reference velocity also had to be reduced. The vehicle could drive at 90mph with 0ms latency, but couldn't go above 75 without driving off the road. Refactoring my code to reduce slow copying of data from one vector to another could also reduce the latency of the entire system if I were to spend more time on it. All messages are supressed to reduce latency.

## Difficulty

It was difficult going into this project with no introduction to the libaries we were using for solving MPC. The fg_eval class and MPC class are very tightly coupled. I spent a lot of time redesigning the system so that the functionality of one class was separated from the other, and spent just as much time putting it back. The two objects end up sharing a lot of member variables in order to avoid a lot of global variables, but I'm still not happy with how it's designed. The fg_eval ()operator and its use in solve is very unintuitive and caused me to waste a lot of time trying to figure out what was going on.

## Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```sh
    git clone https://github.com/uWebSockets/uWebSockets
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.

* **Ipopt and CppAD:** Please refer to [this document](https://github.com/udacity/CarND-MPC-Project/blob/master/install_Ipopt_CppAD.md) for installation instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.
