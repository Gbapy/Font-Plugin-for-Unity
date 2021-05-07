#include "stdafx.h"
#include "ImgCvt.h"
#include "primitive.h"
#include "libFontArt.h"

HBITMAP generateBitmap(std::string text, char *curFont, bool bBold, bool bUnderline, bool bItalic) {
	vector<string> segs;
	int mxLine = 0;
	int lnCount = 1;
	int lnLength = 0;
	int st = 0;
	for (UINT i = 0; i < text.size(); i++) {
		if (text.at(i) == '\n') {
			lnLength = 0;
			lnCount++;
			segs.push_back(text.substr(st, i));
			st = i + 1;
		}
		lnLength++;
		if (lnLength > mxLine) mxLine = lnLength;
	}
	segs.push_back(text.substr(st));
	int width = 300 * mxLine;
	int height = 300 * lnCount;

	HDC		hdc = GetDC(NULL);
	HGDIOBJ previousSelectedHandle = NULL;
	HFONT	hFont = NULL;
	HFONT	previousFont = NULL;
	HBITMAP cropped = NULL;
	HDC		compatibleDeviceContext = NULL;
	HBITMAP bitmapHandle = NULL;
	RECT	rect = RECT();

	
	compatibleDeviceContext = CreateCompatibleDC(hdc);
	bitmapHandle = CreateCompatibleBitmap(hdc, width, height);

	if (bitmapHandle == NULL) goto _final;
	previousSelectedHandle = SelectObject(compatibleDeviceContext, bitmapHandle);

	hFont = CreateFont(300, 0, 0, 0, bBold ? 900 : 300, bItalic, bUnderline, FALSE, DEFAULT_CHARSET, OUT_OUTLINE_PRECIS,
		CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH, TEXT(curFont));

	if (hFont == NULL) goto _final;

	SetTextColor(compatibleDeviceContext, RGB(255, 255, 255));
	SetBkMode(compatibleDeviceContext, TRANSPARENT);
	previousFont = (HFONT)SelectObject(compatibleDeviceContext, hFont);
	for (UINT i = 0; i < segs.size(); i++) {
		rect.left = 0; rect.right = width * (int)(segs[i].size()); rect.top = i * 300; rect.bottom = (i + 1) * 300;
		DrawText(compatibleDeviceContext, TEXT(segs[i].c_str()), (int)segs[i].size(), &rect, DT_EDITCONTROL);
	}
	//TextOut(compatibleDeviceContext, 0, 0, TEXT(text.c_str()), text.size());

	SelectObject(compatibleDeviceContext, previousFont);
	SelectObject(compatibleDeviceContext, previousSelectedHandle);

_final:
	if (compatibleDeviceContext) DeleteDC(compatibleDeviceContext);
	if (hFont) DeleteObject(hFont);
	if (hdc) DeleteDC(hdc);
	return bitmapHandle;
}

vector<BYTE> ToPixels(HBITMAP hBitmap, int &width, int &height)
{
	BITMAP Bmp = { 0 };
	BITMAPINFO Info = { 0 };
	std::vector<BYTE> Pixels = std::vector<BYTE>();
	HDC DC = CreateCompatibleDC(NULL);
	HBITMAP OldBitmap = (HBITMAP)SelectObject(DC, hBitmap);

	std::memset(&Info, 0, sizeof(BITMAPINFO)); //not necessary really..

	GetObject(hBitmap, sizeof(Bmp), &Bmp);

	Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	Info.bmiHeader.biWidth = width = Bmp.bmWidth;
	Info.bmiHeader.biHeight = height = Bmp.bmHeight;
	Info.bmiHeader.biPlanes = 1;
	Info.bmiHeader.biBitCount = Bmp.bmBitsPixel;
	Info.bmiHeader.biCompression = BI_RGB;
	Info.bmiHeader.biSizeImage = ((width * Bmp.bmBitsPixel + 31) / 32) * 4 * height;

	Pixels.resize(Info.bmiHeader.biSizeImage);
	GetDIBits(DC, hBitmap, 0, height, &Pixels[0], &Info, DIB_RGB_COLORS);
	SelectObject(DC, OldBitmap);

	height = std::abs(height);
	DeleteDC(DC);
	return Pixels;
}

bool crop_extend(UCHAR *src, int width, int height,
	UCHAR **dst, int *newWidth, int *newHeight, int margin) {
	int mnI = height;
	int mnJ = width;
	int mxI = 0;
	int mxJ = 0;

	for (int i = 0; i < height; i++) {
		UCHAR *ptr = src + i * width;
		for (int j = 0; j < width; j++) {
			if (ptr[j] == 0) continue;
			if (i < mnI) mnI = i;
			if (j < mnJ) mnJ = j;
			if (i > mxI) mxI = i;
			if (j > mxJ) mxJ = j;
		}
	}

	*newWidth = mxJ - mnJ + 1 + margin * 2;
	*newHeight = mxI - mnI + 1 + margin * 2;

	if (*newWidth < 1 || *newHeight < 1) return false;
	*dst = (UCHAR *)malloc((*newWidth) * (*newHeight));
	if (*dst == NULL) return false;
	memset(*dst, 0, (*newWidth) * (*newHeight));

	for (int i = mnI; i <= mxI; i++) {
		UCHAR *pSrc = src + i * width;
		UCHAR *pDst = *dst + (i + margin - mnI) * (*newWidth);
		for (int j = mnJ; j <= mxJ; j++) {
			pDst[j - mnJ + margin] = pSrc[j];
		}
	}
	return true;
}

void extractContourPixels(UCHAR *src, int width, int height) {
	UCHAR *tmp = (UCHAR *)malloc(width * height);
	memset(tmp, 0, width * height);

	for (int i = 1; i < height - 1; i++) {
		UCHAR *p1 = src + (i - 1) * width;
		UCHAR *p2 = src + i * width;
		UCHAR *p3 = src + (i + 1) * width;
		UCHAR *ptrDst = tmp + i * width;

		for (int j = 1; j < width - 1; j++) {
			if (p2[j] == 0) continue;
			int n = 0;

			n += p1[j] == 0 ? 1 : 0;
			n += p2[j - 1] == 0 ? 1 : 0;
			n += p2[j + 1] == 0 ? 1 : 0;
			n += p3[j] == 0 ? 1 : 0;
			if (n > 0) ptrDst[j] = 255;
		}
	}

	memcpy(src, tmp, width * height);
	free(tmp);
}

int getClosedContour(UCHAR *src, int width, vector<Vec2i> *path, int pi, int pj, bool apartFlag, int nCount) {
	if (nCount > 1000) return -1;
	int ci = -1; int cj = -1;
	for (int k = 0; k < 8; k++) {
		int i1 = pi + neighbor[k].val[0];
		int j1 = pj + neighbor[k].val[1];
		if (src[i1 * width + j1] == 0) continue;
		vector<Vec2i> curPath;
		bool flag = abs(ci - path->at(0).val[0]) > 1 || abs(cj - path->at(0).val[1]) > 1 ?
			true : false;
		ci = i1; cj = j1;
		path->push_back(Vec2i(ci, cj));
		src[ci * width + cj] = 0;
		if (apartFlag && abs(ci - path->at(0).val[0]) < 2 && abs(cj - path->at(0).val[1]) < 2) {
			return 1;
		}
		int ret = getClosedContour(src, width, path, ci, cj, flag, nCount + 1);
		if (ret != 0) return ret;
		if (ret == 0) path->pop_back();
	}
	return 0;
}

TERMINAL twistPoint(double x, double y, TERMINAL tCenter, int width, double xTwist) {
	if (xTwist == 0) return TERMINAL(x, y);
	double radius = tCenter.y - y;
	double length = x - tCenter.x;
	double delta = length * (xTwist * 0.5f) / ((double)width / 2.0f);
	return TERMINAL(tCenter.x + radius * sin(delta * PI / 180.0f), tCenter.y - radius * cos(delta * PI / 180.0f));
}

void twistPath(vector<Vec2i> *path, TERMINAL tCenter, int width, double xTwist) {
	if (xTwist == 0) return;
	for (int i = 0; i < (int)(path->size()); i++) {
		double x = path->at(i).val[1];
		double y = path->at(i).val[0];
		TERMINAL t = twistPoint(x, y, tCenter, width, xTwist);
		path->at(i).val[1] = (int)(t.x);
		path->at(i).val[0] = (int)(t.y);
	}
}

void twistNode(vector<SUBNODE> *nodes, double yTwist) {
	if (yTwist == 0.0f) return;

	double mnX = M_INFINITE;
	double mxX = -M_INFINITE;
	double mnY = M_INFINITE;
	double mxY = -M_INFINITE;
	for (UINT i = 0; i < nodes->size(); i++) {
		for (UINT j = 0; j < nodes->at(i).vertex.size(); j++) {
			VERTEX v = nodes->at(i).vertex[j];
			if (mnX > v.val[0]) mnX = v.val[0];
			if (mxX < v.val[0]) mxX = v.val[0];
			if (mnY > v.val[1]) mnY = v.val[1];
			if (mxY < v.val[1]) mxY = v.val[1];
		}
	}

	int segCount = (int)(fabs(yTwist) / 5.0f) + 1;
	double radius = (mxX - mnX) / (yTwist * PI / 180.0f);
	double stp = (mxX - mnX) / (double)segCount;
	for (int n = 1; n < segCount; n++) {
		double x = mnX + stp * n;
		for (UINT i = 0; i < nodes->size(); i++) {
			int m = (int)(nodes->at(i).index.size());
			for (int j = 0; j < m; j++) {
				TRIANGLE tr = nodes->at(i).index[j];
				VERTEX v1 = nodes->at(i).vertex[tr.idx[0]];
				VERTEX v2 = nodes->at(i).vertex[tr.idx[1]];
				VERTEX v3 = nodes->at(i).vertex[tr.idx[2]];
				v2 = doMinus(v2, v1);
				v3 = doMinus(v3, v1);
				VERTEX v0 = doCross(v2, v3);
				doNormalize(&v0);
				vector<VERTEX> sharePoints;
				int unsharedIndex = -1;
				for (int k = 0; k < 3; k++) {
					int p = k + 1 == 3 ? 0 : k + 1;
					v1 = nodes->at(i).vertex[tr.idx[k]];
					v2 = nodes->at(i).vertex[tr.idx[p]];
					if (isConflict(v1, v2, x, &v3) == 1) {
						sharePoints.push_back(v3);
					}
					else {
						unsharedIndex = k;
					}
				}
				if (sharePoints.size() == 2 && unsharedIndex != -1) {
					int p = unsharedIndex + 2 >= 3 ? unsharedIndex - 1 : unsharedIndex + 2;
					int k = unsharedIndex + 1 == 3 ? 0 : unsharedIndex + 1;
					v1 = sharePoints[0];
					v2 = sharePoints[1];
					v3 = nodes->at(i).vertex[tr.idx[p]];
					v1 = doMinus(v1, v3);
					v2 = doMinus(v2, v3);
					v1 = doCross(v1, v2);
					doNormalize(&v1);
					v1 = doMinus(v1, v0);
					if (doMagnitude(v1) < EP) {
						nodes->at(i).vertex.push_back(sharePoints[0]);
						nodes->at(i).vertex.push_back(sharePoints[1]);
					}
					else {
						nodes->at(i).vertex.push_back(sharePoints[1]);
						nodes->at(i).vertex.push_back(sharePoints[0]);
					}
					int vc = (int)(nodes->at(i).vertex.size());
					nodes->at(i).index[j].idx[0] = vc - 2;
					nodes->at(i).index[j].idx[1] = vc - 1;
					nodes->at(i).index[j].idx[2] = tr.idx[p];
					nodes->at(i).index.push_back(TRIANGLE(vc - 2, tr.idx[k], vc - 1));
					nodes->at(i).index.push_back(TRIANGLE(vc - 2, tr.idx[unsharedIndex], tr.idx[k]));
				}
			}
		}
	}

	double ct = (mnX + mxX) / 2.0f;
	double length = (mxX - mnX) / 2.0f;
	for (UINT i = 0; i < nodes->size(); i++) {
		int m = (int)(nodes->at(i).vertex.size());
		for (int j = 0; j < m; j++) {
			VERTEX v = nodes->at(i).vertex[j];
			v.uv[0] = (v.val[0] - mnX) / (mxX - mnX);
			v.uv[1] = (v.val[1] - mnY) / (mxY - mnY);
			int s1 = (int)((v.val[0] - mnX) / stp);
			double x1 = mnX + s1 * stp;
			double x2 = x1 + stp;
			double rate = (v.val[0] - x1) / stp;
			double r = radius + v.val[2];
			x1 = x1 - ct;
			x2 = x2 - ct;
			x1 = x1 * yTwist * 0.5f / length;
			x2 = x2 * yTwist * 0.5f / length;
			double z1 = radius - r * cos(x1 * PI / 180.0f);
			double z2 = radius - r * cos(x2 * PI / 180.0f);
			x1 = ct + r * sin(x1 * PI / 180.0f);
			x2 = ct + r * sin(x2 * PI / 180.0f);
			v.val[0] = x1 + (x2 - x1) * rate;
			v.val[2] = z1 + (z2 - z1) * rate;
			nodes->at(i).vertex[j] = v;
		}
	}
	for (UINT i = 0; i < nodes->size(); i++) {
		for (UINT j = 0; j < nodes->at(i).index.size(); j++) {
			TRIANGLE tr = nodes->at(i).index[j];
			tr.flip();
			nodes->at(i).index[j] = tr;
		}
	}
}

void generateFlatMeshes(SHAPE *cont, SUBNODE *sb) {
	vector<TERMINAL> verts;
	vector<int>		vertIndexs;
	int prevVertCount = (int)sb->vertex.size();

	for (UINT i = 0; i < cont->prims.size(); i++) {
		verts.push_back(cont->prims[i]->terms[0]);
		VERTEX v = VERTEX(cont->prims[i]->terms[0].x, cont->prims[i]->terms[0].y, 0);

		sb->vertex.push_back(v);
		vertIndexs.push_back(i);
	}

	while (verts.size() > 2) {
		int vCount = (int)vertIndexs.size();

		for (int i = 0; i < (int)vertIndexs.size() - 1; i++) {
			int p = (int)i - 1 < 0 ? (int)vertIndexs.size() - 1 : i - 1;
			int j = i + 1;
			int k = i + 2 >= (int)vertIndexs.size() ? 0 : i + 2;
			int m = i + 3 >= (int)vertIndexs.size() ? i + 3 - (int)vertIndexs.size() : i + 3;
			LINE ln1 = LINE(verts[vertIndexs[i]], verts[vertIndexs[j]]);
			LINE ln2 = LINE(verts[vertIndexs[j]], verts[vertIndexs[k]]);
			LINE ln3 = LINE(verts[vertIndexs[k]], verts[vertIndexs[i]]);
			SHAPE shp;
			TERMINAL t;
			int retIndex;

			shp.prims.push_back(&ln1);
			shp.prims.push_back(&ln2);
			shp.prims.push_back(&ln3);

			if (cont->getIntersection(&ln3, &t, &retIndex)) continue;

			VERTEX v1 = ln1.getPositiveDirection();
			VERTEX v2 = ln3.getNegativeDirection();

			v1 = doCross(v1, v2);
			if (v1.val[2] > 0 && shp.isInsidePoint(verts[m]) == false && shp.isInsidePoint(verts[p]) == false) {
				TRIANGLE tr;
				tr.idx[0] = prevVertCount + vertIndexs[k];
				tr.idx[1] = prevVertCount + vertIndexs[i];
				tr.idx[2] = prevVertCount + vertIndexs[j];
				sb->index.push_back(tr);
				vertIndexs.erase(vertIndexs.begin() + j);
			}
		}
		if (vCount < 3) break;
		if (vCount == (int)(vertIndexs.size())) break;
	}
	verts.clear();
	vertIndexs.clear();
}

void unifyContour(vector<SHAPE> *contGroup) {
	vector<SHAPE> shps;
	TERMINAL t;
	int retIdx;
	for (UINT n = 0; n < contGroup->size(); n++) {
		contGroup->at(n).isValid = true;
	}
	for (UINT n = 1; n < contGroup->size(); n++) {
		if (contGroup->at(n).isValid == false) continue;
		int primCount = (int)(contGroup->at(n).prims.size());
		for (int i = 0; i < primCount; i++) {
			TERMINAL t1 = contGroup->at(n).prims[i]->terms[0];
			if (t1.isValid == false) continue;
			bool flag = true;

			for (UINT m = 0; m < contGroup->size(); m++) {
				if (contGroup->at(m).isValid == false) continue;
				if (m == n) continue;
				for (UINT j = 0; j < contGroup->at(m).prims.size(); j++) {
					TERMINAL t2 = contGroup->at(m).prims[j]->terms[0];
					if (t2.isValid == false) continue;
					LINE ln = LINE(t1, t2);
					flag = true;
					for (UINT p = 0; p < contGroup->size(); p++) {
						if (contGroup->at(p).isValid == false) continue;
						if (contGroup->at(p).getIntersection(&ln, &t, &retIdx) == true) {
							flag = false;
							break;
						}
					}
					if (flag) {
						contGroup->at(n).prims[i]->terms[0].isValid = false;
						contGroup->at(m).prims[j]->terms[0].isValid = false;
						t1.isValid = false;
						t2.isValid = false;
						contGroup->at(m).prims.insert(contGroup->at(m).prims.begin() + j,
							new LINE(t1, t2));
						int pCount = (int)(contGroup->at(n).prims.size());
						SHAPE shp;
						for (int p = 0; p < pCount; p++) {
							int k = i + p >= pCount ? i + p - pCount : i + p;
							shp.prims.push_back(contGroup->at(n).prims[k]);
						}
						contGroup->at(m).prims.insert(contGroup->at(m).prims.begin() + j,
							shp.prims.begin(), shp.prims.end());
						contGroup->at(m).prims.insert(contGroup->at(m).prims.begin() + j,
							new LINE(t2, t1));
						contGroup->at(n).isValid = false;
						break;
					}
				}
				if (flag) break;
			}
			if (flag) break;
		}
	}
}

void smoothenContour(vector<SHAPE> *conts) {
	//Smoothen contours
	for (int n = 0; n < (int)conts->size(); n++) {
		for (int i = 0; i < (int)(conts->at(n).prims.size()); i++) {
			int j = i + 1 >= (int)(conts->at(n).prims.size()) ? 0 : i + 1;
			double m1 = conts->at(n).prims[i]->getMagnitude();
			double m2 = conts->at(n).prims[j]->getMagnitude();
			VERTEX g1 = conts->at(n).prims[i]->getGradient();
			VERTEX g2 = conts->at(n).prims[j]->getGradient();
			VERTEX d1 = conts->at(n).prims[i]->getPositiveDirection();
			VERTEX d2 = conts->at(n).prims[j]->getPositiveDirection();
			d2 = doCross(d1, d2);
			if (d2.val[2] < 0) {
				g1.inverseDirection(); g2.inverseDirection();
			}
			double m = m1 < m2 ? m1 : m2;
			m = m / 3.0f;
			TERMINAL t1 = conts->at(n).prims[i]->getNegativeOffset(m / m1);
			TERMINAL t2 = conts->at(n).prims[j]->getPositiveOffset(m / m2);
			conts->at(n).prims[i]->terms[1] = t1;
			conts->at(n).prims[j]->terms[0] = t2;
			conts->at(n).prims.insert(conts->at(n).prims.begin() + j, new LINE(t1, t2));
			i++;
		}
	}
}

void buildBackWall(double depth, vector<SUBNODE> *nodes) {
	int nCount = (int)(nodes->size());
	for (int i = 0; i < nCount; i++) {
		SUBNODE sb;
		for (UINT j = 0; j < nodes->at(i).vertex.size(); j++) {
			VERTEX v = nodes->at(i).vertex[j];
			v.val[2] = depth;
			sb.vertex.push_back(v);
		}
		for (UINT j = 0; j < nodes->at(i).index.size(); j++) {
			TRIANGLE tr = nodes->at(i).index[j];
			int t = tr.idx[1]; tr.idx[1] = tr.idx[2]; tr.idx[2] = t;
			sb.index.push_back(tr);
		}
		nodes->push_back(sb);
	}
}

void buildFrontWall(vector<SHAPE> *conts, vector<SUBNODE> *nodes) {
	vector<vector<SHAPE>> contourGroup;
	//Regrouping the contours
	for (UINT pi = 0; pi < conts->size(); pi++) {
		vector<SHAPE> contGroup;
		if (conts->at(pi).isValid == false) continue;
		contGroup.push_back(conts->at(pi));
		for (UINT pj = pi + 1; pj < conts->size(); pj++) {
			if (conts->at(pi).isInsideShape(&(conts->at(pj))) == true) {
				bool isAlreadyContained = false;

				for (int k = 1; k < (int)contGroup.size(); k++) {
					if (contGroup[k].isInsideShape(&(conts->at(pj))) == true) {
						isAlreadyContained = true;
						break;
					}
				}
				if (isAlreadyContained == false) {
					conts->at(pj).isValid = false;
					conts->at(pj).isPositive = false;
					contGroup.push_back(conts->at(pj));
				}
			}
		}
		for (UINT i = 1; i < contGroup.size(); i++) {
			contGroup[i].reversePrimitive();
		}
		contourGroup.push_back(contGroup);
	}

	nodes->clear();
	for (UINT n = 0; n < contourGroup.size(); n++) {
		SUBNODE sb;
		unifyContour(&(contourGroup[n]));
		generateFlatMeshes(&(contourGroup[n][0]), &sb);
		nodes->push_back(sb);
	}
}

void buildSideWall(vector<SHAPE> *conts, double depth, vector<SUBNODE> *nodes) {
	SUBNODE sb;
	for (int i = 0; i < (int)(conts->size()); i++) {
		int p = (int)(conts->at(i).prims.size());
		int vCount = (int)(sb.vertex.size());
		for (int j = 0; j < p; j++) {
			int k = j + 1 >= p ? 0 : j + 1;
			TERMINAL t = conts->at(i).prims[j]->terms[0];
			sb.vertex.push_back(VERTEX(t.x, t.y, conts->at(i).depth));
			sb.vertex.push_back(VERTEX(t.x, t.y, depth));
			TRIANGLE tr1 = TRIANGLE(vCount + j * 2, vCount + j * 2 + 1, vCount + k * 2 + 1);
			TRIANGLE tr2 = TRIANGLE(vCount + j * 2, vCount + k * 2 + 1, vCount + k * 2);
			if (conts->at(i).isPositive == false) {
				tr1.flip(); tr2.flip();
			}
			sb.index.push_back(tr1);
			sb.index.push_back(tr2);
		}
	}
	if (sb.vertex.size() > 0 && sb.index.size() > 0) nodes->push_back(sb);
}

void simplifyContour(vector<SHAPE> *conts) {
	for (UINT n = 0; n < conts->size(); n++) {
		for (int i = 0; i < (int)conts->at(n).prims.size(); i++) {
			int k = i + 1 == (int)conts->at(n).prims.size() ? 0 : i + 1;
			int j = i - 1 >= 0 ? i - 1 : (int)conts->at(n).prims.size() - 1;
			if (conts->at(n).prims[i]->getMagnitude() < 2.0) {
				conts->at(n).prims[j]->terms[1] = conts->at(n).prims[k]->terms[0];
				conts->at(n).prims.erase(conts->at(n).prims.begin() + i);
				i--;
			}
			else {
				VERTEX v1 = conts->at(n).prims[i]->getPositiveDirection();
				VERTEX v2 = conts->at(n).prims[k]->getPositiveDirection();
				v1 = doCross(v1, v2);
				if (fabs(v1.val[2]) < 0.1f) {
					conts->at(n).prims[i]->terms[1] = conts->at(n).prims[k]->terms[1];
					conts->at(n).prims.erase(conts->at(n).prims.begin() + k);
					i--;
				}
			}
		}
		conts->at(n).isValid = true;
	}
}

void turnUpsideDownNode(vector<SUBNODE> *nodes, int height) {
	for (UINT i = 0; i < nodes->size(); i++) {
		for (UINT j = 0; j < nodes->at(i).vertex.size(); j++) {
			VERTEX v = nodes->at(i).vertex[j];
			nodes->at(i).vertex[j].val[1] = (double)height - v.val[1];
		}
	}
}

void PutColor(UCHAR *src, int width, int height, double x, double y, UCHAR r, UCHAR g, UCHAR b) {
	if (x < 0 || x >= width) return;
	if (y < 0 || y >= height) return;
	src[(int)y * width * 3 + (int)x * 3] = r;
	src[(int)y * width * 3 + (int)x * 3 + 1] = g;
	src[(int)y * width * 3 + (int)x * 3 + 2] = b;
}

void DrawLine(UCHAR *src, int width, int height, double x1, double y1, double x2, double y2, UCHAR r, UCHAR g, UCHAR b) {
	double m = sqrtf((float)((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)));
	double stpX = (x2 - x1) / (double)m;
	double stpY = (y2 - y1) / (double)m;
	double sX = x1;
	double sY = y1;
	for (int i = 0; i < m; i++) {
		sX += stpX;
		sY += stpY;
		PutColor(src, width, height, sX, sY, r, g, b);
	}
}

void writeContourImage(vector<SHAPE> *conts, int width, int height) {
	UCHAR *dmp = (UCHAR *)malloc(width * height * 3);
	memset(dmp, 0, width * height * 3);
	for (int n = 0; n < (int)(conts->size()); n++) {
		for (int i = 0; i < (int)(conts->at(n).prims.size()); i++) {
			DrawLine(dmp, width, height, conts->at(n).prims[i]->terms[0].x,
				height - conts->at(n).prims[i]->terms[0].y, conts->at(n).prims[i]->terms[1].x,
				height - conts->at(n).prims[i]->terms[1].y, 0, 0, 255);
		}
	}

	generateBitmapImage(dmp, height, width, "D:\\contour.bmp");
	free(dmp);
}

void writeMonoImage(UCHAR *src, int width, int height) {
	UCHAR *dmp = (UCHAR *)malloc(width * height * 3);
	for (int i = 0; i < height; i++) {
		UCHAR *p1 = dmp + i * width * 3;
		UCHAR *p2 = src + i * width;
		for (int j = 0; j < width; j++) {
			p1[j * 3] = p2[j];
			p1[j * 3 + 1] = p2[j];
			p1[j * 3 + 2] = p2[j];
		}
	}
	generateBitmapImage(dmp, height, width, "D:\\result.bmp");
	free(dmp);
}

void centeringMesh(vector<SUBNODE> *nodes) {
	VERTEX mn = VERTEX(M_INFINITE, M_INFINITE, M_INFINITE);
	VERTEX mx = VERTEX(-M_INFINITE, -M_INFINITE, -M_INFINITE);
	for (UINT i = 0; i < nodes->size(); i++) {
		for (UINT j = 0; j < nodes->at(i).vertex.size(); j++) {
			VERTEX v = nodes->at(i).vertex[j];
			if (v.val[0] > mx.val[0]) mx.val[0] = v.val[0];
			if (v.val[1] > mx.val[1]) mx.val[1] = v.val[1];
			if (v.val[2] > mx.val[2]) mx.val[2] = v.val[2];
			if (v.val[0] < mn.val[0]) mn.val[0] = v.val[0];
			if (v.val[1] < mn.val[1]) mn.val[1] = v.val[1];
			if (v.val[2] < mn.val[2]) mn.val[2] = v.val[2];
		}
	}
	VERTEX ct = VERTEX((mn.val[0] + mx.val[0]) / 2.0f, (mn.val[1] + mx.val[1]) / 2.0f,
		(mn.val[2] + mx.val[2]) / 2.0f);
	for (UINT i = 0; i < nodes->size(); i++) {
		for (UINT j = 0; j < nodes->at(i).vertex.size(); j++) {
			VERTEX v = nodes->at(i).vertex[j];
			v = doMinus(v, ct);
			v.nor = nodes->at(i).vertex[j].nor;
			nodes->at(i).vertex[j] = v;
		}
	}
}

int generateModelFromText(string renderText, char *fontName, bool bBold, bool bUnderline, bool bItalic, float thickness,
	double xTwist, double yTwist) {
	vector<BYTE> pixel;
	int		rawWidth, rawHeight;
	int		width, height;
	UCHAR	*raw = NULL;
	UCHAR	*src = NULL;
	SUBNODE subnode;
	vector<SHAPE> conts;

	HBITMAP hBmp = generateBitmap(renderText, fontName, bBold, bUnderline, bItalic);
	if (hBmp == NULL) 
		return FONTART_FAIL_GENERATE_BMP;
	pixel = ToPixels(hBmp, rawWidth, rawHeight);
	if (rawWidth <= 0 || rawHeight <= 0) 
		return FONTART_FAIL_GENERATE_BMP;

	raw = (UCHAR *)malloc(rawHeight * rawWidth);
	if (raw == NULL) 
		return FONTART_FAIL_ALLOC_MEMORY;

	for (int i = 0; i < rawHeight; i++) {
		UCHAR *p = raw + i * rawWidth;
		for (int j = 0; j < rawWidth; j++) {
			p[j] = pixel[i * rawWidth * 4 + j * 4];
		}
	}
	pixel.clear();
	crop_extend(raw, rawWidth, rawHeight, &src, &width, &height, MARGIN);
	if (src == NULL || width <= 0 || height <= 0) 
		return FONTART_FAIL_CROP;

	extractContourPixels(src, width, height);
	writeMonoImage(src, width, height);

	double length = (double)(width - MARGIN * 2);
	double radius = xTwist != 0 ? length / (xTwist * PI / 180.0f) : 0;
	TERMINAL tCenter;
	if (xTwist > 0) {
		tCenter = TERMINAL((double)width / 2.0f, (double)height - MARGIN + radius);
	}
	else if (xTwist < 0) {
		tCenter = TERMINAL((double)width / 2.0f, (double)MARGIN + radius);
	}

	for (int i = 0; i < height; i++) {
		UCHAR *p = src + i * width;
		for (int j = 0; j < width; j++) {
			if (p[j] == 0) continue;
			int ci = i; int cj = j;
			int primCount = 0;

			SHAPE shp;
			vector<Vec2i> path;
			path.push_back(Vec2i(ci, cj));
			p[j] = 0;

		_retry:
			int ret = getClosedContour(src, width, &path, ci, cj, false, 0);
			if (ret == -1) {
				ci = path[path.size() - 1].val[0];
				cj = path[path.size() - 1].val[1];

				goto _retry;
			}
			twistPath(&path, tCenter, width - MARGIN * 2, xTwist);
			int n = 0;
			int m = (int)path.size();
			for (int pi = 1; pi < m; pi++) {
				VERTEX v_xy = VERTEX(path[pi].val[0] - path[n].val[0], path[pi].val[1] - path[n].val[1], 0);
				double length_xy = GetMagnitude(v_xy);
				bool f = false;

				for (int pj = n; pj < pi; pj++) {
					VERTEX v1_xy = VERTEX(path[pj].val[0] - path[n].val[0], path[pj].val[1] - path[n].val[1], 0);
					VERTEX v2_xy = doCross(v_xy, v1_xy);

					double s = fabs(v2_xy.val[2]) / length_xy;
					if (s > 1.5f) {
						f = true;
						break;
					}
				}
				if (f) {
					LINE *ln = new LINE(TERMINAL(path[n].val[1], height - path[n].val[0]),
						TERMINAL(path[pi - 1].val[1], height - path[pi - 1].val[0]));
					shp.prims.push_back(ln);
					primCount++;
					n = pi - 1;
				}
			}

			LINE *ln = new LINE(TERMINAL(path[n].val[1], height - path[n].val[0]),
				TERMINAL(path[m - 1].val[1], height - path[m - 1].val[0]));
			shp.prims.push_back(ln);
			primCount++;

			if (primCount > 2) {
				if (shp.prims[primCount - 1]->terms[1].isEqual(shp.prims[0]->terms[0]) == false) {
					shp.prims[primCount - 1]->terms[1] = shp.prims[0]->terms[0];
				}
				shp.isGemhole = false;
				shp.depth = 0;
				shp.makePositive();
				conts.push_back(shp);
			}
		}
	}

	//simplifyContour(&conts);
	//smoothenContour(&conts);
	writeContourImage(&conts, width, height);
	buildFrontWall(&conts, &nodes);
	buildBackWall(thickness * 100.0f, &nodes);
	buildSideWall(&conts, thickness * 100.0f, &nodes);
	twistNode(&nodes, yTwist);
	doSetVertexNormal(&nodes);
	turnUpsideDownNode(&nodes, height);
	DoScale(&nodes, 1.0f / 100.0f);
	centeringMesh(&nodes);
	for (UINT i = 0; i < conts.size(); i++) {
		conts[i].clear();
	}

	free(raw);
	free(src);
	DeleteObject(hBmp);
	return FONTART_SUCCESS;
}