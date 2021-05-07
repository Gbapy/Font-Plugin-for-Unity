// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the LIBFONTART_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// LIBFONTART_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

#include "ImgCvt.h"

#ifdef LIBFONTART_EXPORTS
#define LIBFONTART_API __declspec(dllexport)
#else
#define LIBFONTART_API __declspec(dllimport)
#endif


// This class is exported from the libFontArt.dll
class LIBFONTART_API ClibFontArt {
public:
	ClibFontArt(void);
	// TODO: add your methods here.
};

typedef struct _FONT_ {
	char name[255];
}FONT, *PFONT;

enum FONTART_ERROR {
	FONTART_SUCCESS = 0,
	FONTART_FAIL_ALLOC_MEMORY = 1,
	FONTART_FAIL_GENERATE_BMP = 2,
	FONTART_FAIL_CROP = 3
};

extern LIBFONTART_API int nlibFontArt;
extern vector<SUBNODE> nodes;
LIBFONTART_API int fnlibFontArt(void);
int generateModelFromText(string renderText, char *fontName, bool bBold, bool bUnderline, bool bItalic,
	float thickness, double xTwist, double yTwist);