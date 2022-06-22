# P3D Assignment 1 - My RayTracer Group 01

## Group members:
### 93604 - Patrícia Vilão
### 93695 - Catarina Sousa
### 93735 - Maria Ribeiro
<br />

# How to Compile and Link:

We developed the project in Visual Studio 2019, so make sure you have this program installed!

1 -> Click in the file: "Dependencies.exe" (this will extract the necessary files).

2 -> Click in the file: "MyRayTracer.sln" and open it with Visual Studio 2019 when asked.

3 -> Select "x64" on the platform

4 -> Click on "Local Windows Debugger" (the play button)

5 -> Enjoy!

<br />

# Flags:

1 -> Soft Shadows (line 91)

2 -> Anti Aliasing (line 92)

3 -> Schlick Approximation (line 93)

4 -> Depth of Field (line 94)

<br />

### Soft Shadows and Antialiasing:
You can change the value of samples per pixel and the number of lights (lines 34 and 35), but make sure that these values are the same!

<br />

# Extra Work:
1 -> Fuzzy Reflections (line 95)

2 -> Motion Blur (line 96)

<br />

## All these flags start at False, if you want to try them, turn them to True in the corresponding line!


<br />

### Fuzzy Reflections:
If you want to change the roughness value (starts at 0.2) you can change it in line 98.

<br />

### Motion Blur:
You can change the value of the interval of the motion by changing the variable "time_f" in line 36.