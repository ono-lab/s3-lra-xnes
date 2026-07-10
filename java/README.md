# S3-LRA-xNES in Java

## Prerequisites

- Java 17
- Maven

### Version Information

```bash
$ java --version
openjdk 17.0.15 2025-04-15
OpenJDK Runtime Environment (build 17.0.15+6-Ubuntu-0ubuntu122.04)
OpenJDK 64-Bit Server VM (build 17.0.15+6-Ubuntu-0ubuntu122.04, mixed mode, sharing)
```

```bash
$ mvn --version
Apache Maven 3.6.3
Maven home: /usr/share/maven
Java version: 17.0.15, vendor: Ubuntu, runtime: /usr/lib/jvm/java-17-openjdk-amd64
Default locale: *********, platform encoding: UTF-8
OS name: "linux", version: "6.8.0-65-generic", arch: "amd64", family: "unix"
```

## Source Code of S3-LRA-xNES

`src/main/java/sthree_lra_xnes/S3LRAxNES.java`

## Sample Optimization Code with S3-LRA-xNES

`src/main/java/main/Main.java`

## Usage

### Executing the Sample Code

```bash
mvn package
java -cp target/S3LRA-xNES-1.0-SNAPSHOT.jar main.Main # executing the sample optimization code with S3-LRA-xNES
```

### Standard Output Example

Executing `src/main/java/main/Main.java` with the command above will produce an output similar to this:

```text
Success:true
NoOfEvals:9408
Best:9.928525751495464E-9
```

- `Success`: A flag that indicates whether the optimization succeeded (`true`) or failed (`false`).
- `NoOfEvals`: The number of function evaluations performed until the stopping criteria were met.
- `Best`: The best objective value found during the optimization.

### How to Modify the Sample Code

To change the experiment settings (other than the objective function), modify the following section in `src/main/java/main/Main.java`:

```java
        /** The problem dimension. **/
        int dimension = 20;
        /** The initial mean vector. */
        double meanElem = 3;
        /** The initial step-size. */
        double sigma = 2;
        /** The random seed. */
        long seed = 0l;
        /** The evaluation budget */
        long maxEval = dimension == 80 ? (long) 2e7 : (long) 1e7;
        /** The target function value. */
        double criterion = 1e-8;
```

For example, if you change `meanElem=8` and `sigma=7`, the initial distribution is changed to $\mathcal{N}(\bm{8},7^2\bm{I})$.

To use a different objective function, first, implement it as a method that takes a `TCMatrix` object (the solution vector) and returns a `double` (the evaluation value):

```java
    static double YOUR_OBJECTIVE_FUNCTION(TCMatrix x) {
        // Implement your function here
        return 0.0;
    }
```

Then, modify the line where the objective function is called inside the `while` loop:

```java
                    // Change the objective function if you want:
                    double value = SphereFunction(s.getVector());
                    // double value = Rastrigin(s.getVector()); // like this,
                    // double value = Ellipsoid(s.getVector()); // or like this.
```

... like this:

```java
                    // Change the objective function if you want:
                    double value = YOUR_OBJECTIVE_FUNCTION(s.getVector());
```
