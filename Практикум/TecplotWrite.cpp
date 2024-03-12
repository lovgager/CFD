#if(defined(_WIN32))
#include <direct.h>
#endif


void Cabaret_LeftTrans_NoSonic::WriteToTecplot(int step, double time)
{
	stringstream ss;
	struct stat s;
	string DirRes = "out";
	if (stat(DirRes.c_str(), &s) != 0)
	{
#if(defined(_WIN32))
		_mkdir(DirRes.c_str());
#else
		mkdir(DirRes.c_str(), 0700);
#endif
	}

	ss << step;
#if(defined(_WIN32))
	string filename = DirRes + "\\out_" + ss.str() + ".dat";
#else
	string filename = DirRes + "/out_" + ss.str() + ".dat";
#endif  

	ofstream out(filename.c_str());
	if (out.is_open())
	{
		out << "TITLE=\"OUT_CONSVT\"" << endl;
		out << "VARIABLES = \"X\", \"uc\"" << endl;
		out << "ZONE" << endl;
		out << "T = CabLeftNoSonic_FluxCorr_uc, I = " << N << endl;
		out << "ZONETYPE = ORDERED, DATAPACKING = POINT, C = BLUE" << endl;
		out << "STRANDID = 1" << endl;
		out << "SOLUTIONTIME = " << time << endl;
		out << endl;
		for (int i = 0; i < N; i++)
			out << xc[i] << " " << u_consvt[i] << endl;
		out << endl;
		out << "TITLE=\"OUT_FLUX\"" << endl;
		out << "VARIABLES = \"X\", \"uf\"" << endl;
		out << "ZONE" << endl;
		out << "T = CabLeftNoSonic_FluxCorr_uf, I = " << N + 1 << endl;
		out << "ZONETYPE = ORDERED, DATAPACKING = POINT, C = BLUE" << endl;
		out << "STRANDID = 2" << endl;
		out << "SOLUTIONTIME = " << time << endl;
		out << endl;
		for (int i = 0; i <= N; i++)
			out << x[i] << " " << u_flux[i] << endl;
		out.close();
	}
	else
	{
		cerr << "WARNING: THE RESULTS ON TIME " << time << " WERE NOT SAVED, STEP = " << step << ", METHOD = CabLeftNoSonic" << endl;
	}
}