// libFontArt.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "libFontArt.h"
#include <assimp/cexport.h>
#include <assimp/scene.h>
#include <assimp/postprocess.h>


vector<FONT> fonts;
vector<SUBNODE> nodes;
// This is an example of an exported variable
LIBFONTART_API int nlibFontArt=0;

// This is an example of an exported function.
LIBFONTART_API int fnlibFontArt(void)
{
	return 42;
}

// This is the constructor of a class that has been exported.
// see libFontArt.h for the class definition
ClibFontArt::ClibFontArt()
{
	return;
}

void ToSBCS(LPWSTR lpszText)
{
	int j = (int)wcslen(lpszText);
	if (j == 0)
	{
		strcpy((char *)lpszText, "");
		return;
	}
	else
	{
		char *lpszNewText = (char *)malloc(j + 1);
		j = WideCharToMultiByte(CP_ACP, 0L, lpszText, -1, lpszNewText, j + 1, NULL, NULL);
		if (j > 0)
			strcpy((char *)lpszText, lpszNewText);
		else
			strcpy((char *)lpszText, "");
		free(lpszNewText);
	}
}

LIBFONTART_API int initOperation() {
	neighbor[0] = Vec2i(-1, 0);
	neighbor[1] = Vec2i(-1, 1);
	neighbor[2] = Vec2i(0, 1);
	neighbor[3] = Vec2i(1, 1);
	neighbor[4] = Vec2i(1, 0);
	neighbor[5] = Vec2i(1, -1);
	neighbor[6] = Vec2i(0, -1);
	neighbor[7] = Vec2i(-1, -1);

	IDWriteFontCollection* pFontCollection = NULL;
	UINT32 familyCount = 0;

	IDWriteFactory *pDWriteFactory;
	HRESULT hr = DWriteCreateFactory(
		DWRITE_FACTORY_TYPE_SHARED,
		__uuidof(IDWriteFactory),
		reinterpret_cast<IUnknown**>(&pDWriteFactory));
	// Get the system font collection.
	if (!SUCCEEDED(hr))
	{
		return -1;// "Error: DWriteCreateFactory", "Error", MB_ICONHAND | MB_OK);
	}

	hr = pDWriteFactory->GetSystemFontCollection(&pFontCollection);
	if (!SUCCEEDED(hr))
	{
		return -2;// "Error: GetSystemFontCollection", "Error", MB_ICONHAND | MB_OK);
	}
	fonts.clear();
	familyCount = pFontCollection->GetFontFamilyCount();
	for (UINT32 i = 0; i < familyCount; ++i) {
		IDWriteFontFamily* pFontFamily = NULL;
		IDWriteLocalizedStrings* pFamilyNames = NULL;
		UINT32 index = 0;
		BOOL exists = false;
		wchar_t localeName[LOCALE_NAME_MAX_LENGTH];
		wchar_t name[LOCALE_NAME_MAX_LENGTH];
		UINT32 length = 0;

		hr = pFontCollection->GetFontFamily(i, &pFontFamily);
		if (!SUCCEEDED(hr))
		{
			continue;
		}
		hr = pFontFamily->GetFamilyNames(&pFamilyNames);
		if (!SUCCEEDED(hr))
		{
			continue;
		}

		int defaultLocaleSuccess = GetUserDefaultLocaleName(localeName, LOCALE_NAME_MAX_LENGTH);

		// If the default locale is returned, find that locale name, otherwise use "en-us".
		if (defaultLocaleSuccess)
		{
			hr = pFamilyNames->FindLocaleName(localeName, &index, &exists);
		}
		if (SUCCEEDED(hr) && !exists) // if the above find did not find a match, retry with US English
		{
			hr = pFamilyNames->FindLocaleName(L"en-us", &index, &exists);
		}
		if (!exists) index = 0;
		if (!SUCCEEDED(hr))
		{
			continue;
		}
		hr = pFamilyNames->GetStringLength(index, &length);

		if (!SUCCEEDED(hr) || name == NULL || length == 0) {
			continue;
		}
		hr = pFamilyNames->GetString(index, name, length + 1);

		ToSBCS(name);
		FONT font;
		memset(font.name, 0, 255);
		strcpy(font.name, (char *)name);
		fonts.push_back(font);
	}
	return (int)(fonts.size());
}

LIBFONTART_API int getFontList(char *buffer) {
	for (UINT i = 0; i < fonts.size(); i++) {
		char *tmp = (char *)buffer + i * 255;
		tmp[0] = (char)strlen(fonts[i].name);
		tmp++;
		strcpy(tmp, (char *)fonts[i].name);
	}
	return 0;
}

LIBFONTART_API int generateModel(char *text, int txtLength, char *fontName, int fntLength, bool bBold, bool bUnderline, bool bItalic, float thickness,
	double xTwist, double yTwist, int *bufferInfo) {
	char txt[256];
	char fnt[256];
	if (text == NULL || fontName == NULL) return -1;

	memset(txt, 0, 256); memset(fnt, 0, 256);
	memcpy(txt, text, txtLength);
	memcpy(fnt, fontName, fntLength);

	for (UINT i = 0; i < nodes.size(); i++) {
		nodes[i].vertex.clear();
		nodes[i].index.clear();
	}
	int ret = generateModelFromText(string(txt), fnt, bBold, bUnderline, bItalic, thickness, xTwist, yTwist);
	bufferInfo[0] = (int)nodes.size();

	for (UINT i = 0; i < nodes.size(); i++) {
		bufferInfo[1] += (int)(nodes[i].vertex.size());
		bufferInfo[2] += (int)(nodes[i].index.size());
	}
	return ret;
}

LIBFONTART_API void getMeshData(int *vBufInfo, int *mBufInfo, float *vBuf, int *mBuf, float *nBuf) {
	int vCount = 0;
	int mCount = 0;
	for (UINT i = 0; i < nodes.size(); i++) {
		vBufInfo[i] = (int)(nodes[i].vertex.size());
		mBufInfo[i] = (int)(nodes[i].index.size());
		for (UINT j = 0; j < nodes[i].vertex.size(); j++) {
			vBuf[vCount * 3] = (float)(nodes[i].vertex[j].val[0]);
			vBuf[vCount * 3 + 1] = (float)(nodes[i].vertex[j].val[1]);
			vBuf[vCount * 3 + 2] = (float)(nodes[i].vertex[j].val[2]);
			nBuf[vCount * 3] = (float)(nodes[i].vertex[j].nor.val[0]);
			nBuf[vCount * 3 + 1] = (float)(nodes[i].vertex[j].nor.val[1]);
			nBuf[vCount * 3 + 2] = (float)(nodes[i].vertex[j].nor.val[2]);
			vCount++;
		}
		for (UINT j = 0; j < nodes[i].index.size(); j++) {
			mBuf[mCount * 3] = nodes[i].index[j].idx[0];
			mBuf[mCount * 3 + 1] = nodes[i].index[j].idx[1];
			mBuf[mCount * 3 + 2] = nodes[i].index[j].idx[2];
			mCount++;
		}
	}
}

int exportObj(vector<SUBNODE> *model, char *fileName, char *ext) {
	if (model->size() == 0) return -2;
	aiScene *scene = new aiScene();

	scene->mRootNode = new aiNode();

	scene->mMaterials = new aiMaterial*[model->size()];
	scene->mNumMaterials = (int)model->size();


	scene->mMeshes = new aiMesh*[model->size()];
	scene->mNumMeshes = (int)model->size();

	scene->mRootNode->mMeshes = new UINT[model->size()];
	scene->mRootNode->mNumMeshes = (int)model->size();

	for (UINT i = 0; i < model->size(); i++) {
		scene->mMaterials[i] = new aiMaterial();

		scene->mMeshes[i] = new aiMesh();
		scene->mMeshes[i]->mMaterialIndex = i;

		auto pMesh = scene->mMeshes[i];
		scene->mRootNode->mMeshes[i] = i;

		pMesh->mVertices = new aiVector3D[model->at(i).vertex.size()];
		pMesh->mNormals = new aiVector3D[model->at(i).vertex.size()];
		pMesh->mNumVertices = (int)model->at(i).vertex.size();

		pMesh->mTextureCoords[0] = new aiVector3D[model->at(i).vertex.size()];
		pMesh->mNumUVComponents[0] = (int)model->at(i).vertex.size();

		for (UINT j = 0; j < model->at(i).vertex.size(); j++)
		{
			VERTEX v = model->at(i).vertex[j];

			pMesh->mVertices[j] = aiVector3D((float)v.val[0], (float)v.val[1], (float)v.val[2]);
			pMesh->mTextureCoords[0][j] = aiVector3D((float)v.uv[0], (float)v.uv[1], 1.0f);
			pMesh->mNormals[j] = aiVector3D((float)v.nor.val[0], (float)v.nor.val[1], (float)v.nor.val[2]);
		}

		pMesh->mFaces = new aiFace[model->at(i).index.size()];
		pMesh->mNumFaces = (UINT)(model->at(i).index.size());

		for (UINT j = 0; j < model->at(i).index.size(); j++)
		{
			aiFace &face = pMesh->mFaces[j];
			TRIANGLE t = model->at(i).index[j];
			face.mIndices = new UINT[3];
			face.mNumIndices = 3;

			face.mIndices[0] = t.idx[0];
			face.mIndices[1] = t.idx[1];
			face.mIndices[2] = t.idx[2];
		}
	}
	aiReturn ret = aiExportScene(scene, ext, fileName, 0);

	for (UINT i = 0; i < model->size(); i++) {
		auto pMesh = scene->mMeshes[i];
		free((void *)pMesh->mVertices);
		free((void *)pMesh->mNormals);
		free((void *)pMesh->mTextureCoords[0]);
		free((void *)scene->mMaterials[i]);
		//free((void *)pMesh->mFaces);
		free((void *)pMesh);
	}
	free((void *)scene->mRootNode->mMeshes);
	free((void *)scene->mRootNode);
	free((void *)scene->mMaterials);
	free((void *)scene);
	return ret == aiReturn_SUCCESS ? 0 : -1;
}

LIBFONTART_API int exportMeshData(float *vertBuffer, float *norBuffer,
	int *meshBuffer, int vertCount, int meshCount, char *path, int pathLength) {

	if (vertCount == 0 || meshCount == 0 || vertBuffer == NULL ||
		norBuffer == NULL || meshBuffer == NULL) return -1;
	vector<SUBNODE> mesh;
	VERTEX mx, mn;
	char *pathBuf = (char *)malloc(pathLength + 1);

	if (pathBuf == NULL) return 1;
	memset(pathBuf, 0, pathLength + 1);
	strcpy(pathBuf, path);

	SUBNODE sb;
	vector<SUBNODE> subNode;
	int sCount = -1;

	for (int i = 0; i < vertCount; i++) {
		VERTEX v = VERTEX(vertBuffer[0], vertBuffer[1], -vertBuffer[2]);
		v.nor.val[0] = norBuffer[0]; v.nor.val[1] = norBuffer[1]; v.nor.val[2] = norBuffer[2];
		sb.vertex.push_back(v);
		vertBuffer += 3; norBuffer += 3;
	}

	for (int i = 0; i < meshCount; i++) {
		TRIANGLE tr = TRIANGLE(meshBuffer[0], meshBuffer[1], meshBuffer[2]);
		tr.flip();
		sb.index.push_back(tr);
		meshBuffer += 3;
	}

	mesh.push_back(sb);

	return exportObj(&mesh, pathBuf, "obj");
}