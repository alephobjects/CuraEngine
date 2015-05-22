/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include "infill.h"
#include <clipper/clipper.hpp>

// sqrt(2) helpers for truncated octahedron calculations
#define SQRT2MUL(x) ((30547*(x))/21600)
#define OCTSLEN(x) ((43200*(x))/73747)
#define OCTDLEN(x) ((21600*(x))/30547)

namespace cura {

void generateConcentricInfill(Polygons outline, Polygons& result, int inset_value)
{
    while(outline.size() > 0)
    {
        for (unsigned int polyNr = 0; polyNr < outline.size(); polyNr++)
        {
            PolygonRef r = outline[polyNr];
            result.add(r);
        }
        outline = outline.offset(-inset_value);
    }
}

void generateAutomaticInfill(const Polygons& in_outline, Polygons& result,
                             int extrusionWidth, int lineSpacing,
                             int infillOverlap, double rotation, int posZ)
{
    if (lineSpacing > extrusionWidth * 6)
    {
        generateTroctInfill(in_outline, result, extrusionWidth, lineSpacing,
                           infillOverlap, rotation, posZ);
    }
    if (lineSpacing > extrusionWidth * 4)
    {
        generateGridInfill(in_outline, result, extrusionWidth, lineSpacing,
                           infillOverlap, rotation);
    }
    else
    {
        generateLineInfill(in_outline, result, extrusionWidth, lineSpacing,
                           infillOverlap, rotation);
    }
}

void generateGridInfill(const Polygons& in_outline, Polygons& result,
                        int extrusionWidth, int lineSpacing, int infillOverlap,
                        double rotation)
{
    generateLineInfill(in_outline, result, extrusionWidth, lineSpacing * 2,
                       infillOverlap, rotation);
    generateLineInfill(in_outline, result, extrusionWidth, lineSpacing * 2,
                       infillOverlap, rotation + 90);
}

int compare_int64_t(const void* a, const void* b)
{
    int64_t n = (*(int64_t*)a) - (*(int64_t*)b);
    if (n < 0) return -1;
    if (n > 0) return 1;
    return 0;
}

void generateLineInfill(const Polygons& in_outline, Polygons& result, int extrusionWidth, int lineSpacing, int infillOverlap, double rotation)
{
    Polygons outline = in_outline.offset(extrusionWidth * infillOverlap / 100);
    PointMatrix matrix(rotation);
    
    outline.applyMatrix(matrix);
    
    AABB boundary(outline);
    
    boundary.min.X = ((boundary.min.X / lineSpacing) - 1) * lineSpacing;
    int lineCount = (boundary.max.X - boundary.min.X + (lineSpacing - 1)) / lineSpacing;
    vector<vector<int64_t> > cutList;
    for(int n=0; n<lineCount; n++)
        cutList.push_back(vector<int64_t>());

    for(unsigned int polyNr=0; polyNr < outline.size(); polyNr++)
    {
        Point p1 = outline[polyNr][outline[polyNr].size()-1];
        for(unsigned int i=0; i < outline[polyNr].size(); i++)
        {
            Point p0 = outline[polyNr][i];
            int idx0 = (p0.X - boundary.min.X) / lineSpacing;
            int idx1 = (p1.X - boundary.min.X) / lineSpacing;
            int64_t xMin = p0.X, xMax = p1.X;
            if (p0.X > p1.X) { xMin = p1.X; xMax = p0.X; }
            if (idx0 > idx1) { int tmp = idx0; idx0 = idx1; idx1 = tmp; }
            for(int idx = idx0; idx<=idx1; idx++)
            {
                int x = (idx * lineSpacing) + boundary.min.X + lineSpacing / 2;
                if (x < xMin) continue;
                if (x >= xMax) continue;
                int y = p0.Y + (p1.Y - p0.Y) * (x - p0.X) / (p1.X - p0.X);
                cutList[idx].push_back(y);
            }
            p1 = p0;
        }
    }
    
    int idx = 0;
    for(int64_t x = boundary.min.X + lineSpacing / 2; x < boundary.max.X; x += lineSpacing)
    {
        qsort(cutList[idx].data(), cutList[idx].size(), sizeof(int64_t), compare_int64_t);
        for(unsigned int i = 0; i + 1 < cutList[idx].size(); i+=2)
        {
            if (cutList[idx][i+1] - cutList[idx][i] < extrusionWidth / 5)
                continue;
            PolygonRef p = result.newPoly();
            p.add(matrix.unapply(Point(x, cutList[idx][i])));
            p.add(matrix.unapply(Point(x, cutList[idx][i+1])));
        }
        idx += 1;
    }
}

  void generateTroctInfill(const Polygons& in_outline, Polygons& result, int extrusionWidth, int lineSpacing, int infillOverlap, double rotation, int posZ)
  {
    Polygons outline = in_outline.offset(extrusionWidth * infillOverlap / 100);
    PointMatrix matrix(rotation);
    outline.applyMatrix(matrix);
    AABB boundary(outline);

    // ignore infill for areas smaller than line spacing
    if((abs(boundary.min.X - boundary.max.X) + abs(boundary.min.Y - boundary.max.Y)) < lineSpacing){
      return;
    }

    // fix to normalise against diagonal infill
    lineSpacing = lineSpacing * 2;

    uint64_t Zscale = SQRT2MUL(lineSpacing);

    int offset = abs(posZ % (Zscale) - (Zscale/2)) - (Zscale/4);
    boundary.min.X = ((boundary.min.X / lineSpacing) - 1) * lineSpacing;
    boundary.min.Y = ((boundary.min.Y / lineSpacing) - 1) * lineSpacing;

    unsigned int lineCountX = (boundary.max.X - boundary.min.X + (lineSpacing - 1)) / lineSpacing;
    unsigned int lineCountY = (boundary.max.Y - boundary.min.Y + (lineSpacing - 1)) / lineSpacing;
    int rtMod = int(rotation / 90) % 2;
    // with an odd number of lines, sides need to be swapped around
    if(rtMod == 1){
      rtMod = (lineCountX + int(rotation / 90)) % 2;
    }

    // draw non-horizontal walls of octohedrons
    Polygons po;
    PolygonRef p = po.newPoly();
    for(unsigned int ly=0; ly < lineCountY;){
      for(size_t it = 0; it < 2; ly++, it++){
        int side = (2*((ly + it + rtMod) % 2) - 1);
        int y = (ly * lineSpacing) + boundary.min.Y + lineSpacing / 2 - (offset/2 * side);
        int x = boundary.min.X-(offset/2);
        if(it == 1){
          x = (lineCountX * (lineSpacing)) + boundary.min.X + lineSpacing / 2 - (offset/2);
        }
        p.add(Point(x,y));
        for(unsigned int lx=0; lx < lineCountX; lx++){
          if(it == 1){
            side = (2*((lx + ly + it + rtMod + lineCountX) % 2) - 1);
            y = (ly * lineSpacing) + boundary.min.Y + lineSpacing / 2 + (offset/2 * side);
            x = ((lineCountX - lx - 1) * lineSpacing) + boundary.min.X + lineSpacing / 2;
            p.add(Point(x+lineSpacing-abs(offset/2), y));
            p.add(Point(x+abs(offset/2), y));
          } else {
            side = (2*((lx + ly + it + rtMod) % 2) - 1);
            y = (ly * lineSpacing) + boundary.min.Y + lineSpacing / 2 + (offset/2 * side);
            x = (lx * lineSpacing) + boundary.min.X + lineSpacing / 2;
            p.add(Point(x+abs(offset/2), y));
            p.add(Point(x+lineSpacing-abs(offset/2), y));
          }
        }
        x = (lineCountX * lineSpacing) + boundary.min.X + lineSpacing / 2 - (offset/2);
        if(it == 1){
          x = boundary.min.X-(offset/2);
        }
        y = (ly * lineSpacing) + boundary.min.Y + lineSpacing / 2 - (offset/2 * side);
        p.add(Point(x,y));
      }
    }
    // Generate tops / bottoms of octohedrons
    if(abs((abs(offset) - Zscale/4)) < (extrusionWidth/2)){
      uint64_t startLine = (offset < 0) ? 0 : 1;
      uint64_t coverWidth = OCTSLEN(lineSpacing);
      vector<Point> points;
      for(size_t xi = 0; xi < (lineCountX+1); xi++){
        for(size_t yi = 0; yi < (lineCountY); yi += 2){
          points.push_back(Point(boundary.min.X + OCTDLEN(lineSpacing)
                                 + (xi - startLine + rtMod) * lineSpacing,
                                 boundary.min.Y + OCTDLEN(lineSpacing)
                                 + (yi + (xi%2)) * lineSpacing
                                 + extrusionWidth/2));
        }
      }
      uint64_t order = 0;
      for(Point pp : points){
        PolygonRef p = po.newPoly();
        for(size_t yi = 0; yi <= coverWidth; yi += extrusionWidth) {
          if(order == 0){
            p.add(Point(pp.X, pp.Y + yi));
            p.add(Point(pp.X + coverWidth + extrusionWidth, pp.Y + yi));
          } else {
            p.add(Point(pp.X + coverWidth + extrusionWidth, pp.Y + yi));
            p.add(Point(pp.X, pp.Y + yi));
          }
          order = (order + 1) % 2;
        }
      }
    }
    // intersect with outline polygon(s)
    Polygons pi = po.intersection(outline);
    // Hack to add intersection to result. There doesn't seem
    // to be a direct way to do this
    for(unsigned int polyNr=0; polyNr < pi.size(); polyNr++) {
      PolygonRef p = result.newPoly(); //  = result.newPoly()
      for(unsigned int i=0; i < pi[polyNr].size(); i++) {
        Point p0 = pi[polyNr][i];
        p.add(matrix.unapply(Point(p0.X,p0.Y)));
      }
    }
}

}//namespace cura
