/*//EWave Implementation by Soumitra Goswami
Author: Soumitra Goswami
Date: 5/08/2018
Last Edited: 5/10/2018

Description: Original paper by Dr.Jerry Tessendorf. This is my implementation of it. The program intends to generate
			 a simulation of an interactive wave propagation as a 2D displacement map. 
*/
#include <Windows.h>
#include <omp.h>
#include <GL\glew.h>
#include <GL\freeglut.h>
#include <glm/glm.hpp>
#include <iostream>
#include <string>
#include <ctime>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include "eWaveSim.h"

using namespace std;

int iWidth, iHeight, iChannels;
float *pixmap;
float *imageFile;
unsigned char *outputFile;

int frame;
float brightnessScale;

int paint_mode, display_mode;
enum { PAINT_OBSTRUCTION, PAINT_SOURCE, PAINT_DIVERGENCE, PAINT_COLOR };
enum { OBSTRUCTION_DISPLAY, HEIGHT_DISPLAY, VELOCITY_DISPLAY, PRESSURE_DISPLAY, DIVERGENCE_DISPLAY };

#define SWAP(x,y) {float *temp=x;x=y;y=temp; }
#define BRUSH_SIZE 31
float source_brush[BRUSH_SIZE][BRUSH_SIZE];
float obstruction_brush[BRUSH_SIZE][BRUSH_SIZE];
int nXGrids, nYGrids;
float xLength, yLength;
float dXSize, dYSize;

int x_mouse_prev, y_mouse_prev;

float timeStep = 0.03f;
float scaling_factor = 1.0f;
float *source_height, *source_obstruction;


eWaveSim *sim;

//INPUTS
float gravity[2] = { 0,40.0f };

bool pause_sim = false;
bool capture_screen = false;
//-------------------------------IO FUNCTIONS------------------------------------//

//************************************
// Method:    readImage
// FullName:  readImage
// Access:    public 
// Returns:   void
// Parameter: const char * fileName - Name of the file to be read .
// Parameter: float * & map - Pixmap to store the image.

// Description: Reads in a image from file and stores it in a pixmap.
//************************************
void readImage(const char* fileName, float *&map)
{
	int xres, yres, channels;
	float *inImage = stbi_loadf(fileName, &xres, &yres, &channels, 0);
	if (!inImage) { return; }
	iWidth = xres; // Assigns to global width parameter
	iHeight = yres; // Assigns to global height parameter
	iChannels = channels; // Assigns to global channel parameter
	map = new float[xres*yres*channels];
	long index = 0;
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			for (int c = 0; c < channels; c++)
			{
				long iIndex = (i + xres*(yres - j - 1))*channels + c;
				//cout << inImage[index] << endl;
				(map)[iIndex] = inImage[index++];
				//cout << map[iIndex] << endl;
			}
		}
	}
	//delete[] inImage;
}

//************************************
// Method:    writeImage
// FullName:  writeImage
// Access:    public 
// Returns:   void
// Parameter: const char * fileName - output file name
// Parameter: float * & map - Pixmap used to write the file. In this program we use jpegs

// Description: Uses STB Method to write the file to location of the program
//************************************

void writeImage(const char* fileName, float *&map)
{
	//JPEG quality for the written image
	int quality = 100;
	for (int j = 0; j < iHeight; j++)
	{
		for (int i = 0; i < iWidth; i++)
		{
			int index = i + j*(iWidth);
			int out_index = (iHeight - j - 1) * iWidth + i;
			//Converts colors from floating point to unsigned ints
			int r = int(map[index * 3 + 0] * 255);
			int g = int(map[index * 3 + 1] * 255);
			int b = int(map[index * 3 + 2] * 255);
			
			// Making sure the pixmap is in range from 0-255
			if (r < 0) r = 0; if (r > 255) r = 255;
			if (g < 0) g = 0; if (g > 255) g = 255;
			if (b < 0) b = 0; if (b > 255) b = 255;

			outputFile[out_index * 3 + 0] = r;
			outputFile[out_index * 3 + 1] = g;
			outputFile[out_index * 3 + 2] = b;
		}
	}

	stbi_write_jpg(fileName, iWidth, iHeight, 3, outputFile, quality);


}

//--------------------------------------------INITIALIZE FUNCTIONS---------------------------------------------------//

//************************************
// Method:    initMaps
// FullName:  Initialize Maps
// Access:    public 
// Returns:   void
// Parameter: float(* & map) - The array to initialize
// Parameter: int size - total size of the array. (Usually = Width*Height)
// Parameter: float value - the value to assign to the entire array
// 
// Description: This function initializes the maps to a floating point value. Used to quickly initialize Arrays.
//************************************

static void initMaps(float(*&map), int size, float value)
{

#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		(map)[i] = value;
	}
}



//************************************
// Method:    initializeBrush
// FullName:  initializeBrush
// Access:    public 
// Returns:   void
// 
// Description: Initializes the painting brush for source and obstruction.
//			    Radius of the brush is determined by the BRUSH_SIZE
//			    Brush uses a Circular Kernal
//************************************

void initializeBrush()
{
	int brush_rad = (BRUSH_SIZE - 1) / 2;
	for (int j = -brush_rad; j <= brush_rad; j++)
	{
		int jj = j + brush_rad;
		float jRatio = (float(brush_rad) - fabs(j)) / float(brush_rad);
		for (int i = -brush_rad; i <= brush_rad; i++)
		{
			int ii = i + brush_rad;
			float iRatio = (float(brush_rad) - fabs(i)) / float(brush_rad);
			float radius2 = (jRatio*jRatio + iRatio*iRatio) / 2.0;
			//Generating a Circular Kernal.
			source_brush[ii][jj] = pow(radius2, 0.5);
			//Circular Kernal with 0 at the center.
			obstruction_brush[ii][jj] = 1.0 - pow(radius2, 1.0 / 4.0);
		}
	}
}

//************************************
// Method:    initializeFields
// FullName:  initializeFields
// Access:    public 
// Returns:   void
// 
// Description: Quick Function to initialize all the required maps to be used later for the program.
//************************************

void initializeFields()
{
	int size = iWidth*iHeight;
	initMaps(source_height, size, 0.0);
	initMaps(source_obstruction, size, 1.0);
}


//************************************
// Method:    printMap
// FullName:  printMap
// Access:    public 
// Returns:   void

// Parameter: float * & map - The pixmap to be printed out
// Parameter: int x - Width of the pixmap
// Parameter: int y - Height of the pixmap
// Parameter: int c - No of channels in the pixmap
// Parameter: char * name - Name of the pixmap
// 
// Description: Quick Debugging tool to print out the contents of a pixmap.
//************************************
void printMap(float *&map, int x, int y, int c, char* name)
{
	for (int j = 0; j < y; j++)
	{
		for (int i = 0; i < x; i++)
		{
			for (int k = 0; k < c; k++)
			{
				long index = (i + x*j)*c + k;
				cout << name << " Pixel(" << i << "," << j << "): " << map[index] << endl;
			}
		}
	}
}

//************************************
// Method:    setNbCores
// FullName:  setNbCores
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: int nb - Number of cores to set
// 
// Description: A quick function to set the number of threads for parallel processing
//************************************
void setNbCores(int nb)
{
	omp_set_num_threads(nb);
}

// -------------------------------------IMAGE AND PAINT FUNCTIONS------------------------------------//

//************************************
// Method:    displayImage
// FullName:  displayImage
// Access:    public 
// Returns:   void

//Description: Displaying the pixmap on the screen. Display switches between Obstruction Display 
//			   and non Obstruction display based on the user's mode selection.
//************************************

void displayImage(void)
{
	if (frame != 1)
	{
		//printMap(source_obstruction, iWidth, iHeight, 1, "obsMap");
		//printMap(old_r, iWidth, iHeight, 1);
	}
	for (int j = 0; j < iHeight; j++)
	{
		#pragma omp parallel for
		for (int i = 0; i < iWidth; i++)
		{
			int index = (i + iWidth*j) * 3;
			float r, g, b;
			r = g = b = 0;
			if (display_mode == OBSTRUCTION_DISPLAY)
			{
				float col = 0.5 * ((sim->height(i,j) / scaling_factor) + 1.0)*source_obstruction[index/3];
				r = col;
				g = col;
				b = col;
			}
			else if (display_mode == HEIGHT_DISPLAY)
			{
				float col = 0.5 * ((sim->height(i, j) / scaling_factor) + 1.0);
				r = col;
				g = col;//sim.densityField[index / 3];
				b = col;//sim.densityField[index / 3];
			}
			else if (display_mode == VELOCITY_DISPLAY)//NOT IN USE
			{
				int i = (index / 3) * 2;
				r = 0;// fabs(sim.velField[i]);
				g = 0;//fabs(sim.velField[i + 1]);
				b = 0.0;
			}
			else if (display_mode == PRESSURE_DISPLAY)// NOT IN USE
			{
				r = 0;
				g = 0;
				b = 0;
			}
			
			pixmap[index + 0] = r * brightnessScale;
			pixmap[index + 1] = g * brightnessScale;
			pixmap[index + 2] = b * brightnessScale;
			//cout << "Pixel: " << (index / 3) << " r: " << r << " g: " << g << " b: " << b << endl;
		}
	}
}

//************************************
// Method:    scaleBrightness
// FullName:  scaleBrightness
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: float value - The value to use to scale brightness.
// 
// Description: Quick method to scale the brightness on the screen.
//************************************

void scaleBrightness(float value)
{
	brightnessScale *= value;
	cout << "BRIGHTNESS: " << brightnessScale << endl;
}

//************************************
// Method:    paintScreen
// FullName:  paintScreen
// Access:    public 
// Returns:   void
// 
// Parameter: int x - x position on screen
// Parameter: int y - y position on screen
// 
// Description: This function adds the height and obstruction sources based on the user's mouse inputs.  
//************************************

void paintScreen(int x, int y)
{
	int brush_radius = (BRUSH_SIZE - 1) / 2;
	int xstart = x - brush_radius;
	int ystart = y - brush_radius;

	if (xstart < 0) { xstart = 0; }
	if (ystart < 0) { ystart = 0; }

	int xend = x + brush_radius;
	int yend = y + brush_radius;
	if (xend >= iWidth) { xend = iWidth - 1; }
	if (yend >= iHeight) { yend = iHeight - 1; }

	if (paint_mode == PAINT_SOURCE)
	{
		for (int ix = xstart; ix <= xend; ix++)
		{
			for (int iy = ystart; iy <= yend; iy++)
			{
				int index = ix + iWidth*(iHeight - iy - 1);
				imageFile[3 * index + 0] = sim->height(index);
				imageFile[3 * index + 1] = sim->height(index);
				imageFile[3 * index + 2] = sim->height(index);
				source_height[index] += source_brush[ix - xstart][iy - ystart];

			}
		}
	}
	if (paint_mode == PAINT_OBSTRUCTION)
	{
		for (int ix = xstart; ix <= xend; ix++)
		{
			for (int iy = ystart; iy <= yend; iy++)
			{
				int index = ix + iWidth*(iHeight - iy - 1);
				imageFile[3 * index + 0] *= obstruction_brush[ix - xstart][iy - ystart];
				imageFile[3 * index + 1] *= obstruction_brush[ix - xstart][iy - ystart];
				imageFile[3 * index + 2] *= obstruction_brush[ix - xstart][iy - ystart];
				source_obstruction[index] *= glm::clamp(obstruction_brush[ix - xstart][iy - ystart], 0.0f, 1.0f);
			}
		}
	}
	return;
}


// ---------------------------GLUT FUNCTIONS ---------------------------------------------------------//

void gRender(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(iWidth, iHeight, GL_RGB, GL_FLOAT, pixmap);
	glutSwapBuffers();
}


void gIdleState(void)
{
	clock_t current_ticks = clock();
	displayImage();

	if (!pause_sim) {
		//paintScreen(256, 400);
		sim->propogate(source_height, source_obstruction,timeStep);
		
		initMaps(source_height, (iWidth*iHeight), 0);
	}
	glutPostRedisplay();
	clock_t delta_ticks = clock() - current_ticks;
	if (capture_screen)
	{
		string advection;
		string dispframe = to_string(frame);
		if (frame < 1000) { dispframe = "0" + dispframe; }
		if (frame < 100) { dispframe = "0" + dispframe; }
		if (frame < 10) { dispframe = "0" + dispframe; }
		string fName = "Renders/SG_eWave_" + advection + "_" + dispframe + ".jpg";

		writeImage(fName.c_str(), pixmap);
		cout << "Writing Frame: " << dispframe << endl;
	}
	
	frame++;
	if (frame % 24 == 0)
	{
		float FPS = 0;
		if (delta_ticks > 0)
			FPS = CLOCKS_PER_SEC / delta_ticks;
		cout << "fps: " << FPS << endl;
	}
}


void gKeyboardControls(unsigned char key, int x, int y)

{
	switch (key)
	{
	case '-': case '_':
		scaleBrightness(0.9f);
		break;

	case '+': case '=':
		scaleBrightness(1.0f / 0.9f);
		break;
	case '[': case '{':
		timeStep += 0.05f;
		cout << "timeStep: " << timeStep << endl;
		break;
	case ']': case '}':
		timeStep -= 0.05f;
		cout << "timeStep: " << timeStep << endl;
		break;
	case ',': case '<':
		gravity[1] *= 0.5;
		cout << "gravity: " << gravity[1] << endl;
		break;
	case '.': case '>':
		gravity[1] *= 2;
		cout << "gravity: " << gravity[1] << endl;
		break;
	case 'c': case 'C':
		if (capture_screen == true) {
			capture_screen = false;
			cout << "CAPTURE OFF" << endl;
		}
		else {
			capture_screen = true;
			cout << "CAPTURE ON" << endl;
		}
		break;
	case 'r':
		brightnessScale = 1.0;
		break;
	case 'm': case 'M':
		display_mode++;
		display_mode = display_mode % 2;
		if (display_mode == OBSTRUCTION_DISPLAY) cout << "DISPLAY_OBSTRUCTION" << endl;
		else if (display_mode == HEIGHT_DISPLAY) cout << "DON'T DISPLAY OBSTRUCTION" << endl;
		break;
	case ' ':
		if (pause_sim)
		{
			cout << "SIMULATION STARTED" << endl;
			pause_sim = false;
		}
		else
		{
			cout << "SIMULATION PAUSED" << endl;
			pause_sim = true;
		}
		break;
	case 's':
		paint_mode = PAINT_SOURCE;
		cout << "PAINTING SOURCE" << endl;
		break;
	case 'o': case 'O':
		paint_mode = PAINT_OBSTRUCTION;
		cout << "PAINTING OBSTRUCTION" << endl;
		break;

	default:
		break;

	}

}

void gMouseDown(int button, int state, int x, int y)
{
	if (button != GLUT_LEFT_BUTTON) { return; }
	if (state != GLUT_DOWN) { return; }
	x_mouse_prev = x;
	y_mouse_prev = y;
	paintScreen(x, y);
}

void gMouseMove(int x, int y)
{
	x_mouse_prev = x;
	y_mouse_prev = y;
	paintScreen(x, y);
}

void PrintUsage()
{
	cout << "FLUID_paint keyboard choices\n";
	cout << "s			turns on painting source strength\n";
	cout << "o			turns on painting obstructions\n";
	cout << "+/-		increase/decrease brightness of display\n";
	cout << "r			resets brightness to default\n";
	cout << "c			toggles screen capture on/off\n";
	cout << "'['or ']'  increase/decrease timestep\n";
	cout << "',' or '.' increase/decrease Gravity\n";
	cout << "'m' or 'M' change debug mode(Color, Density, Velocity, Pressure)\n";
}

int glutFunctions(int argc, char *argv[])
{
	// Initialize GLUT
	glutInit(&argc, argv);
	// Set up some memory buffers for our display
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	// Set the window size
	glutInitWindowSize(iWidth, iHeight);
	// Create the window with the title "SG Fluid"
	glutCreateWindow("SG Fluid");

	glClearColor(1, 1, 1, 1);
	glutDisplayFunc(&gRender);
	glutIdleFunc(&gIdleState);
	glutKeyboardFunc(&gKeyboardControls);
	glutMouseFunc(&gMouseDown);
	glutMotionFunc(&gMouseMove);

	GLenum err = glewInit();
	if (GLEW_OK != err) {
		fprintf(stderr, "GLEW error");
		return 1;
	}

	glutMainLoop();
	return 0;
}

int main(int argc, char* argv[]) {
	frame = 1;
	brightnessScale = 1;
	setNbCores(4);
	iWidth = 512;
	iHeight = 512;
	sim = new eWaveSim();

	//constexpr int blah = sizeof(eWaveSim);
	//readImage("grumpy.jpg", imageFile);
	
	nXGrids = iWidth - 2;
	nYGrids = iHeight - 2;
	sim->initFields(nXGrids + 2, nYGrids + 2);
	xLength = (float)iWidth;
	yLength = (float)iHeight;


	dXSize = xLength / float(nXGrids + 2);
	dYSize = yLength / float(nYGrids + 2);


	int size = iWidth*iHeight;
	int grids = (nXGrids + 2)*(nYGrids + 2);

	pixmap = new float[size * 3];
	cout << "Initializing pixmap" << endl;
	outputFile = new unsigned char[size * 3];
	imageFile = new float[size * 3];

	source_height = new float[size];
	source_obstruction = new float[size];
	initializeFields();
	//printMap(source_obstruction, iWidth, iHeight, 1, "obsMap");

	initializeBrush();
	paint_mode = PAINT_SOURCE;
	display_mode = OBSTRUCTION_DISPLAY;
	PrintUsage();

	int r = glutFunctions(argc, argv);
	delete sim;
	return r;
}