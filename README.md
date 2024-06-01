# CSE305 - Project

## Steps to Set Up and Compile the Project

### Step 1: Extract the ImageMagick Library

```bash
tar -xzvf ImageMagick.tar.gz
```

### Step 2: Build and Set Up the Library

```bash
make
```

At this stage, the library has been built and set up, ready to use.

### Step 3: Find the Include Path

To compile the project, you need to locate the `Magick++.h` header file. Run the following command to find the include path:

```bash
find $HOME/ImageMagick -name "Magick++.h"
```

Copy the result of this command, which is the path to `Magick++.h`.

### Step 4: Compile and Run the Project

Use the path found from the previous step and replace it with `<found_include_path>` in the command below. This command compiles the `main.cpp` file with the necessary include and library paths:

```bash
g++ -std=c++11 -o main main.cpp -I<found_include_path> -L$HOME/ImageMagick/lib -lMagick++-7.Q16HDRI -lMagickCore-7.Q16HDRI -lMagickWand-7.Q16HDRI
```

Finally, execute the compiled program:

```bash
./main <algrithm_number>
```

where `<algrithm_number>` is set to 1 for sequential_simulation, 2 for parallel_step_simulation, 3 for parallel_distinc_simulation, 4 for parallel_combined_simulation steps.
