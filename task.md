
In the two folders above you can find ImageD11 and pyFAI git checkouts. These are both python codes for processing diffraction data. Treat them as 
read only, you can only write in this folder where you find task.md

There is some confusion and ongoing discussion about how to convert calibration parameters between the two codes.

In pyFAI, a small json formatted 'poni' file saves the geometry and the maths should be described well within the project. There are some notions 
for conversion to ImageD11 geometry.

In ImageD11, a 'par' file is used within transform.py. Spatial distortion is applied before any geometry parameters are used.

In ImageD11 image flips are done using o11, o12, o21, o22 matrix elements. All 8 possible image flips are allowed.

In pyFAI then 4 possible image orientations (1,2,3,4) are allowed.

The issue is about mapping a par file to a poni file depending on the tilts.

A correct mapping should give the same image of "twotheta" or "q" values from both programs.

The mapping between azimuth angle (eta in ImageD11) should be made clear.

There is an off-by-0.5 pixel difference between the codes. How this 0.5 varies with image flip is unclear.

To begin with, assume there is no spatial distortion.

Tasks:

0: rewrite this task.md into a clean PLAN.md for an LLM to follow. Ask about anything that is unclear to you.

1: review both repos. locate the geometrical descriptions and also divine them from the code and docs. Identify any contradictions. Code is the 
source of truth

2: work out any potential mathematical mappings to get from a par file to a poni file and back again. Write this out to 'mapping.md' and include 
both algebraic syntax and simple python functions to illustrate the mathematical calculations

3: write a pair of functions to convert a poni file to a par file and a par file to a poni file. These take care of the formats on disc and io.

4: ensure the functions can do a correct round trip in both directions

5: develop clear, simple, test cases that explore all the possible flips allowed. These should use a strongly tilted detector (rx,ry,rz all !=0). They 
should match the tth or q values and make clear the mapping between azimuthal values. 

