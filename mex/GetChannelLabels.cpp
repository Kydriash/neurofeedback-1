#include "lsl_cpp.h"
#include "matrix.h"
#include "mex.h"
#include <iostream>
using namespace std;

/* The gateway function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
if (nlhs !=2)
{
 mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "Two inputs required.");
}
char* predicate = reinterpret_cast<char*>(plhs[0]);
char* value = reinterpret_cast<char*>(plhs[1]);
std::vector<lsl::stream_info> results = lsl::resolve_stream(predicate, value);
	lsl::stream_inlet inlet(results[0]);

	// get the full stream info (including custom meta-data) and dissect it
	lsl::stream_info inf = inlet.info();
	//cout << "The stream's XML meta-data is: " << endl;
	//cout << inf.as_xml() << endl;
	//cout << "The manufacturer is: " << inf.desc().child_value("manufacturer") << endl;
	//cout << "The cap circumference is: " << inf.desc().child("cap").child_value("size") << endl;
	//cout << "The channel labels are as follows:" << endl;
	lsl::xml_element ch = inf.desc().child("channels").child("channel");
	mwSize channel_count = inf.channel_count();
	char** ch_labels = new char*[channel_count];
	for (mwSize k = 0; k < channel_count; k++) {
		//cout << "  " << ch.child_value("label") << endl;
		char* channel = new char[10];
		strcpy(channel, ch.child_value("label"));
		ch_labels[k] = channel;
		//strcpy(ch_labels[k], ch.child_value("label"));
		ch = ch.next_sibling();
	}

	const char** const_ch_labels = (const char**)ch_labels;
        
	mwSize m = channel_count;
	nrhs = channel_count;
	plhs[0] =  mxCreateCharMatrixFromStrings(m, const_ch_labels);

}

