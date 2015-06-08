typdef struct particle
    {
    float mass;
    float xpos, ypos, zpos;
    float xvel, yvel, zvel;
    } particle;

 (*readfile(
for(i = 0; i<num_part; i++)
    {
    const int got = fscanf(f, "%f, %f, %f, %f, %f, %f, %f", inArray[i].mass, inArray[i].xpos, inArray[i].ypos, inArray[i].zpos, inArray[i].xvel, inArray[i].yvel, inArray[i].zvel);
    }
