#include <MeshReconstruction.h>
#include <IO.h>
#include <math.h>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <omp.h>

using namespace MeshReconstruction;
using namespace std;

int main(int argc, char* argv[])
{
	int cubesRes = 50;
	int frameNum = 10;
	int THREAD_NUM = 16;
	int saveObj = 0; int correctTest = 0;

	opterr = 0; int c; int index;

	//read the command-line args
	while ((c = getopt (argc, argv, "con:t:f:")) != -1) {
		switch (c) {
			case 'c':
				break;
			case 'o':
				saveObj = 1;
				break;
			case 'n':
				cubesRes = stoi(optarg);
        			break;
			case 't':
				THREAD_NUM = stoi(optarg);
        			break;
			case 'f':
				frameNum = stoi(optarg);
				break;
      		case '?':
        		if (optopt == 'n' || optopt == 'f')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
				exit(0);
			default:
				exit(0);
		}
	}

	for (index = optind; index < argc; index++)
		printf ("Non-option argument %s\n", argv[index]);

	omp_set_num_threads(THREAD_NUM);
	using std::chrono::high_resolution_clock;
	using std::chrono::time_point;
	using std::chrono::duration;

	printf ("resolution (n) = %d, number of frames (f) = %d, threads = %d\n", cubesRes, frameNum, THREAD_NUM);

	//initialize the sdf of the desired shape
	Fun3s linkSdf = [](Vec3 const& pos)
	{
		double radBig = 0.5;
		double radSmall = 0.2;
		double le = 0.3;
		Vec3 q{ pos.x, max(abs(pos.y) - le, 0.0), pos.z };
		return Norm(Norm(q.x, q.y) - radBig, q.z) - radSmall;
	};

	//initialize the twist transformation function
	FunTr opTwist = [](Fun3s const& sdf, Vec3 const& pos, double twist)
	{
		double c = cos(twist*pos.y);
		double s = sin(twist*pos.y);
		Vec3  q{ c*pos.z - s*pos.x, pos.y, s*pos.z + c*pos.x };
		return sdf(q);
	};

	double maxTwist = 5.0;
	double twist = 0.0;
	Rect3 domain;
	domain.min = { -1.0, -1.0, -1.0 };
	domain.size = { 2.0, 2.0, 2.0 };
	Vec3 cubeSize = domain.size*(1.0/cubesRes);
	int frame;

	time_point<high_resolution_clock> start = high_resolution_clock::now();

	//the main time loop
	#pragma omp parallel for num_thread(THREAD_NUM) firstprivate(twist)
	for (frame = 0; frame < frameNum; frame++)
	{
		int thr_id = omp_get_thread_num();
		for (int k=0; k<thr_id; k++){
			twist += 1.0/double(k) * maxTwist;
		}

		//run the main marching cube function to construct a mesh
		Mesh mesh = MarchCube(linkSdf, opTwist, domain, cubeSize, twist);
		if (saveObj)
		{
			//save the object file if told so
			string filename = "p1p2_link_f" + to_string(frame) + "_n" + to_string(cubesRes) + ".obj";
			WriteObjFile(mesh, filename);
		}

	}

	time_point<high_resolution_clock> end = high_resolution_clock::now();
	duration<double> diff = end-start;

	printf ("Time taken: %f sec \n", diff.count());
}