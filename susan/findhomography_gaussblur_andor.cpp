#include <iostream>
#include <string>
#include <math.h>
#include "/home/susan/micromanager1.3/MMDevice/MMDevice.h"
#include "/home/susan/micromanager1.3/MMDevice/DeviceBase.h"
#include "/home/susan/micromanager1.3/MMDevice/ModuleInterface.h"
#include <sstream>
#include <map>
#include <vector>
#include "/home/susan/micromanager1.3/MMCore/MMCore.h"
#include "/home/susan/micromanager1.3/MMCore/Configuration.h"
#include <cvd/image_io.h>
#include <cvd/glwindow.h>
#include <cvd/gl_helpers.h>
#include <cvd/image_convert.h>
#include <cvd/vision.h>
#include <cstdlib>
#include <gvars3/instances.h>
#include <tag/printf.h>
#include <vector>
#include <TooN/LU.h>
#include <TooN/optimization/downhill_simplex.h>
#include <algorithm>
#include <cvd/vector_image_ref.h>
#include <TooN/SVD.h>
#include <cvd/convolution.h>
#include <cvd/image_convert.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <cstdlib>

//This version of the program finds the appropriate homography for a given camera and SLM combination
//This version has most of the debugging printing and pauses removed, and the writing has been somewhat improved.
using namespace std;
using namespace CVD;
using namespace GVars3;
using namespace TooN;

struct calculate_error
{
	vector<Vector<2> > xy, uv;
//	vector<double> x, y, u, v; //though of course we are in fact using the renormalised values
	double operator()(const Vector<8>& a) const
	{
		Matrix<3> homog;
		homog[0] = makeVector(a[0], a[1], a[2]);
		homog[1] = makeVector(a[3], a[4], a[5]);
		homog[2] = makeVector(a[6], a[7], 1);

		double err_sum=0;
		for(unsigned int i=0; i < uv.size(); i++)
		{
			err_sum += norm_sq(uv[i] - project(homog * makeVector(xy[i][0], xy[i][1], 1)));
		}
		return err_sum;
	}
};

template<class C> double length(const C& v)
{
	return sqrt(v*v);
}


template<class C> double simplex_size(const C& s)
{
	return length(s.get_simplex()[s.get_best()] - s.get_simplex()[s.get_worst()]) / length( s.get_simplex()[s.get_best()]);
}

Matrix<3> homography_matrix(const Vector<8>& v)
{
	Matrix<3> homog;
	homog[0] = v.slice<0,3>();
	homog[1] = v.slice<3, 3>();
  	homog[2][0] = v[6];
	homog[2][1] = v[7];
  	homog[2][2] = 1;
	return homog;
}

typedef pair<Display*, Window> Disp;

//Wrapper around CVD images. This sends the data to the xserver using XPutImage
void draw_image(Disp disp, const Image<byte>& im_)
{
	Display* d = disp.first;
	Window w = disp.second;
	int s = DefaultScreen(d);
	Image<Rgb8> im = convert_image(im_);
	XImage *i;
	
	//This will only work on the right kind of display
	//most are the right kind these days
	i = XCreateImage(d, DefaultVisual(d, s), 24, ZPixmap, 0, (char*)im.data(), im.size().x, im.size().y, 8, im.size().x*4);
	
	XPutImage(d, w, DefaultGC(d, s), i, 0, 0, 0, 0, im.size().x, im.size().y);

	i->data = 0;
	XDestroyImage(i);
}


Disp gimme_a_window(ImageRef size)
{
	Display* 	display;
	Window w;
	
	//Open a connection to the X server
	display = XOpenDisplay(0);

	//Create a window
	w = XCreateSimpleWindow(display, DefaultRootWindow(display), 0, 0, size.x, size.y, 0, 0,0);

	XMapWindow(display, w);

	//Synchronize: flush all requests and wait until they've been processed
	//ie wait until the window appears
	XSync(display,1);
	
	return make_pair(display, w);
}

void move_window(Disp d, int x, int y)
{
	XMoveWindow(d.first, d.second, x, y);
	XSync(d.first, 0);
}


int main()
{
	int i,j, iscale=50, jscale=50, NumHomoPics;
	int MINI=12, MINJ=9;
	int MAXI=18, MAXJ=18;
	int xdisplay=1280, ydisplay=1024;
	double biggestmax=0, smallestmax=65535;
	double cutoff=0.3;
	double StartZ, EndZ;
	vector<ImageRef> maxpos(MAXI*MAXJ);
	vector<Image<byte> > all_images(MAXI*MAXJ);
	ImageRef LCOS_size(1280,1024);
	ImageRef picsize;

	//vector<double> xpos;
	//vector<double> ypos;
	vector<Vector<2> > xypos;
	
	//vector<double> upos;
	//vector<double> vpos;
	vector<Vector<2> > uvpos;
	
	vector<double> brightest(MAXI*MAXJ);

	//set the size of the second monitor 1024x768, 
	ImageRef display_size(xdisplay,ydisplay);

	//provisional values for Ilow and Ihigh.  Ihigh will need to be revised later.
	byte Ilow=50;
	byte Ihigh=100;
	
	CMMCore core;
	core.unloadAllDevices();

	cerr << core.getVersionInfo();

	cout << "point 0" << endl;	
	//load camera 
	core.loadDevice("Andor", "Andor", "Andor");

	//initilise all devices
	core.initializeAllDevices();
	
	//core.setProperty("Andor",  "AD_Converter", "2. 14bit"); 

	//choose camera exposure in ms
        core.setExposure(200);

	
	//set up the GL window for display.
	Disp w = gimme_a_window(display_size);
	//move window to the top left hand corner of the second monitor, -1 ypos so the frame doesn't affect the position	
	move_window(w, 1280,-1);



	//generate the test data set
	//generate the starting images here 
	//0-1000 (or 1200 if changed), steps of 200, 0-750 in steps of 150
	//single bright pixel each time
	for(i=MINI;i<MAXI;i++)
	{
		for(j=MINJ;j<MAXJ;j++)
		{
			cout << "point 1" << endl;	
			
			Image<byte> testimage(LCOS_size);
			testimage.fill(0);
			
			for(int yoff = 0; yoff < 4; yoff++)
				for(int xoff = 0; xoff < 4; xoff++)
					testimage[j*jscale +yoff][i*iscale+xoff] = 255;
			//img_save(testimage, sPrintf("test_images2/test2-%03i-%03j.jpg", i, j));
			
			//display an image
			draw_image(w, testimage);	
			sleep(1);		
	
			//take an image
	        	core.snapImage();
		
			//save image
			unsigned short* image_data = (unsigned short*) core.getImage();
			BasicImage<unsigned short> image(image_data, ImageRef(core.getImageWidth(), core.getImageHeight()));
			Image<byte> debug2(image.size());
			int darkest2 = *min_element(image.begin(), image.end());
			int bright2 = *max_element(image.begin(), image.end());

			for(int r=0; r < image.size().y; r++)
			{	
				for(int c=0;  c < image.size().x; c++)
				{
					debug2[r][c] = (byte)floor(255 * (image[r][c]-darkest2+0.1)/(bright2 -darkest2));
				}
			}
	
			img_save(debug2, sPrintf("test_images2/test-%03i-%03j.jpg", i, j));
			
			Image<double> blurred(image.size());
			Image<double> blurred_big(image.size());
			copy(image.begin(), image.end(), blurred.begin());
			copy(image.begin(), image.end(), blurred_big.begin());
			convolveGaussian(blurred, blurred, 1); 
			convolveGaussian(blurred_big, blurred_big, 20); 
			for(int r=0; r < blurred.size().y; r++)
				for(int c=0;  c < blurred.size().x; c++)
					blurred[r][c]-=blurred_big[r][c];
			
			picsize = image.size(); //this is used later for sizing the renormed position vectors
			
			brightest[i+j*MAXI] = *max_element(blurred.begin(), blurred.end());
			int darkest = *min_element(blurred.begin(), blurred.end());
			int bright = *max_element(blurred.begin(), blurred.end());
			Image<byte> debug(blurred.size());
			for(int r=0; r < blurred.size().y; r++)
			{	
				for(int c=0;  c < blurred.size().x; c++)
				{
					debug[r][c] = (byte)floor(255 * (blurred[r][c]-darkest+0.1)/(bright -darkest));
				}
			}
	
			img_save(debug, sPrintf("test_images/test-%03i-%03j.jpg", i, j));
			all_images[i+j*MAXI] = debug;

			maxpos[i+j*MAXI] = blurred.pos(max_element(blurred.begin(), blurred.end()));	
						
			if(brightest[i+j*MAXI] > biggestmax)
				biggestmax = brightest[i+j*MAXI];
			  
			if(brightest[i+j*MAXI]  < smallestmax)
				smallestmax = brightest[i+j*MAXI];
		}
	}

	cout << "biggestmax  " << biggestmax << endl;
	cout << "smallestmax  " << smallestmax << endl;


	for(i=MINI;i<MAXI;i++)
	{
		for(j=MINJ;j<MAXJ;j++)
		{
			//find the spots in the images taken, remember to control for the fact that the spots 
			//may not be visible in all the images 
			//second condition is for LED, where you are looking at much smaller intensity variations and so have to filter more carefully
			if(brightest[i+j*MAXI] > (biggestmax - biggestmax*cutoff)) 
			//(brightest[i+j*6] > (biggestmax - biggestmax*0.01)&& brightest[i+j*6] < (smallestmax + smallestmax*0.04))
			{	
 				//xpos.push_back((double)i*iscale + 1);			
				//ypos.push_back((double)j*jscale + 1);			
				xypos.push_back(makeVector((double)i*iscale + 1, (double)j*jscale + 1));
				
				cout << "xpos " << (double)i*iscale + 1 << endl;		
				cout << "ypos " << (double)j*jscale + 1 << endl;		
	
				//upos.push_back((double)maxpos[i+j*MAXI].x);
				//vpos.push_back((double)maxpos[i+j*MAXI].y);
				uvpos.push_back(makeVector((double)maxpos[i+j*MAXI].x,(double)maxpos[i+j*MAXI].y));

				Image<Rgb<byte> > tmp = convert_image(all_images[i+j*MAXI]);
				tmp[maxpos[i+j*MAXI]] = Rgb<byte>(0,0,255);
				tmp[maxpos[i+j*MAXI]+ImageRef(1,0)] = Rgb<byte>(0,0,255);
				tmp[maxpos[i+j*MAXI]-ImageRef(1,0)] = Rgb<byte>(0,0,255);
				tmp[maxpos[i+j*MAXI]+ImageRef(0,1)] = Rgb<byte>(0,0,255);
				tmp[maxpos[i+j*MAXI]-ImageRef(0,1)] = Rgb<byte>(0,0,255);
				img_save(tmp, sPrintf("test_images/spot_find-%03i-%03j.jpg", i, j));
				
			}
			else if(brightest[i+j*MAXI] < (biggestmax - biggestmax*cutoff))
				cout << "Not bright enough!" << endl;
			else if(brightest[i+j*MAXI] > (smallestmax + smallestmax*0.04))
				cout << "Too bright!" << endl;
			else
				cout << " WTF !!!!!!!!!!!????????????????????????" << endl;

		}
	}


	//Rescaling section starts here
	//rescaled uv vector = rescale1*(u v 1)T 
	//Do the rescaling by hand rather than setting up matrices since it will be less painful.
	//but set up the rescaling matrices here so it is clear how the parameters are used
	
	//uv rescaling
	double auv = 2.0/((double)picsize.x);
	double buv = -1;
	double cuv = 2.0/((double)picsize.y);
	double duv = -1;

	Matrix<3> rescale1;
	rescale1[0] = makeVector(auv, 0, buv);
	rescale1[1] = makeVector(0, cuv, duv);
	rescale1[2] = makeVector(0, 0, 1);
	
	vector<Vector<2> > uvrenorm;
	//vector<double> urenorm(xypos.size());
	for(i=0;i<xypos.size();i++)
		uvrenorm.push_back(makeVector(auv*uvpos[i][0] + buv,cuv*uvpos[i][1] + duv));
		//urenorm[i] = auv*uvpos[i][0] + buv;

//	vector<double> vrenorm(xypos.size());
//	for(i=0;i<xypos.size();i++)
//		vrenorm[i] = cuv*uvpos[i][1] + duv;

	//xy rescaling
	double axy = 2.0/((double)MAXI*(double)iscale);
	double bxy = -1;
	double cxy = 2.0/((double)MAXJ*(double)jscale);
	double dxy = -1;

	Matrix<3> rescale2;
	rescale2[0] = makeVector(axy, 0, bxy);
	rescale2[1] = makeVector(0, cxy, dxy);
	rescale2[2] = makeVector(0, 0, 1);

	vector<Vector<2> > xyrenorm;
	for(i=0;i<xypos.size();i++)
		xyrenorm.push_back(makeVector(axy*xypos[i][0] + bxy, cxy*xypos[i][1] + dxy));
		//xrenorm[i] = axy*xypos[i][0] + bxy;

//	vector<double> yrenorm(xypos.size());
//	for(i=0;i<xypos.size();i++)
//		yrenorm[i] = cxy*xypos[i][1] + dxy;
	

	//process the test data set to find the homography
	//find an approximate starting point using matrix inversion
	Vector<8> astart;
	Vector<8> start;
	
	//Find an approximate starting point using the Moore-Penrose pseudoinverse
	//note that the previous method only worked if the initial four points were not in a row
	//and also suffered from the problem that x, y, u and v were not renormalised to the range -1,1
	Matrix<Dynamic, 8> atrans(xypos.size()*2, 8); 
	for(i=0;i<xypos.size();i++)
	{
		atrans[i*2]=makeVector(xyrenorm[i][0], xyrenorm[i][1], 1,   0,      0,       0,  -xyrenorm[i][0]*uvrenorm[i][0], -xyrenorm[i][1]*uvrenorm[i][0]);
		atrans[i*2+1]=makeVector(0,     0,         0,  xyrenorm[i][0], xyrenorm[i][1], 1,  -xyrenorm[i][0]*uvrenorm[i][1], -xyrenorm[i][1]*uvrenorm[i][1]);
	}

	//get the pseudo-inverse using singular value decomposition
	SVD<Dynamic, 8> svd(atrans);

	//set up utrans.  Multiply the pseudo inverse of A by uvtrans.
	Vector<> uvrenormtrans(xypos.size()*2); 
	for(i=0;i<xypos.size();i++)
	{
		uvrenormtrans[i*2] = uvrenorm[i][0];
		uvrenormtrans[i*2+1] = uvrenorm[i][1];
	}
	astart = svd.get_pinv() * uvrenormtrans;
	
	Vector<8> a;
	a=astart;
	start=a;
	
	//renormalise to get h from a
	LU<3> lu_decomp(rescale1);
	Matrix<3> inv_rescale1=lu_decomp.get_inverse();
	
	Matrix<3> matrixastart = homography_matrix(astart);
	
	cout << "projected matrix" << endl;
	for(i=0;i<xyrenorm.size(); i++)
		cout << xyrenorm[i][0] << "  " << xyrenorm[i][1] << "          " << uvrenorm[i][0] << "   " << uvrenorm[i][1] << "          " << project(matrixastart * makeVector(xyrenorm[i][0], xyrenorm[i][1], 1)) << endl;
	
	Matrix<3> starth;
	starth=inv_rescale1*matrixastart*rescale2;
	Vector<8> hstart=makeVector(starth[0][0], starth[0][1], starth[0][2], starth[1][0], starth[1][1], starth[1][2],starth[2][0],starth[2][1]);
	Vector<3> xyvector;
	xyvector = makeVector(xypos[0][0], xypos[0][1], 1);	
	Vector<3> uvcalcvector;
	uvcalcvector = starth*xyvector;

	//Use downhill simplex to find the homography translating between the two sets of images.
	//You're optimizing an 8 parameter thing (a homography), so the method needs to take in a Vector<8>
	//The reason for that is that you can now use your class like a function

	calculate_error f;
	
///	f.x = xrenorm; //corrected to take into account normalised position
//	f.y = yrenorm;
//	f.u = urenorm;
//	f.v = vrenorm;
	f.xy = xyrenorm;
	f.uv = uvrenorm;
	
	double err;
	err = f(start);
	cout << "error in program = " << err << endl;

	//you then create a DownhillSimplex class

	DownhillSimplex<8> dh_fixed(f, start, 1);
	
	cout << "simplex size 1 = " << simplex_size(dh_fixed) << endl;

	//simplex_size function terminates when the simplex gets too small
	//Restart simplex to make sure really at minimum?
	while(simplex_size(dh_fixed) > 0.00000000001)
	{
		dh_fixed.iterate(f);
		//cout << "intermediate values " << dh_fixed.get_simplex() << endl;
	}
	
	cout << "simplex size 2 = " << simplex_size(dh_fixed) << endl;
	
	dh_fixed.restart(f, 0.001);

	cout << "Final homography  a" << dh_fixed.get_simplex()[dh_fixed.get_best()] << endl;

	Vector<8> afinal = dh_fixed.get_simplex()[dh_fixed.get_best()];
	

	Matrix<3> matrixa = homography_matrix(afinal);
	
	//convert a to h
	Matrix<3> finalh;
	finalh=inv_rescale1*matrixa*rescale2;

	cout << "projected matrix for hstart" << endl;
	for(i=0;i<xyrenorm.size(); i++)
		cout << xypos[i][0] << "  " << xypos[i][1] << "          " << uvpos[i][0] << "   " << uvpos[i][1] << "          " << project(finalh * makeVector(xypos[i][0], xypos[i][1], 1)) << endl;
	
	cout << "Final homography  h" << endl << finalh <<endl;	
	
}
