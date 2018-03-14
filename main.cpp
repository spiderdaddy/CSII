
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstring>


// Include GLEW
#include <GL/glew.h>
#include <GL/freeglut.h>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;

#include <shader.hpp>
#include <texture.hpp>
#include <controls.hpp>

#define WINDOW_TITLE_PREFIX "Gravity"

struct XY {
    float x, y, z, w;
};

struct RGBA {
    float r, g, b, a;
};

struct SegmentVertices {
    XY v1;
    XY v2;
    XY v3;
    XY v4;
    XY v5;
    XY v6;

};

struct SegmentColours {
    RGBA c1;
    RGBA c2;
    RGBA c3;
    RGBA c4;
    RGBA c5;
    RGBA c6;

};

struct P {
    double v, m;
};

struct Segment {
    double r, theta, x, y;
    double vr, vt;
    double ar, at;
    double area;
    double m;
    int n[8] = { -1, -1, -1, -1, -1, -1, -1, -1}; //nearest neighbours
    std::vector<P> pr;
    std::vector<P> pt;
};

#define G 6.6e-11
int levels = 6;
int num_points = pow(2, levels);
double AU = 1.49e08;
double inner_radius =  0.39 * AU;
double outer_radius = 40.00 * AU;

#define M_EARTH   1e22
#define M_SUN     2e30
#define M_JUPITER 1e27

double stellar_mass = M_SUN;
double escape_mass = 0;

#define STAR
#define SELF

std::vector<SegmentVertices> segmentVertices;
std::vector<SegmentColours> segmentColours;
std::vector<Segment> segment;
std::vector<Segment> newSegment;


int
        CurrentWidth = 1000,
        CurrentHeight = 1000,
        WindowHandle = 0;

unsigned FrameCount = 0;

GLuint  VaoId,
        VboAreaId,
        VboPropertyId;

void Initialize(int, char *[]);

void InitWindow(int, char *[]);

void ResizeFunction(int, int);

void RenderFunction(void);

void TimerFunction(int);

void IdleFunction(void);

void Cleanup(void);

void CreateVBO(void);

void DestroyVBO(void);

void CreateShaders(void);

void DestroyShaders(void);

void PrintPoints();

void MapSegmentToColor();

void CalcSystemMass();

int main(int argc, char *argv[]) {
    Initialize(argc, argv);

    glutMainLoop();

    exit(EXIT_SUCCESS);
}


void Initialize(int argc, char *argv[]) {

    segmentVertices.reserve(num_points*num_points);
    segmentColours.reserve(num_points*num_points);
    segment.reserve(num_points*num_points);
    newSegment.reserve(num_points*num_points);

    srand (time(NULL));

    // populate points
    double theta_step = 2 * M_PI / num_points;
    double theta_step_2 = theta_step / 2;
    double radius_step = ( outer_radius - inner_radius ) / num_points;
    double radius_step_2 = radius_step / 2;

    for (size_t t = 0; t < num_points; t++) {

        double theta = t * theta_step + theta_step_2;

        for (size_t r = 0; r < num_points; r++) {

            double radius = inner_radius + (r * radius_step) + radius_step_2;

            Segment s;
            s.pr.reserve(9);
            s.pt.reserve(9);
            s.r = radius;
            s.theta = theta;
            s.x = radius * cos(theta);
            s.y = radius * sin(theta);
            s.area = 2.0 * radius * radius_step * theta_step;
            s.m = (double)(rand()%50) * s.area; //double(rand() % 100000)/100000.0f * max_segment_mass / 1000.0f;
            s.vr = 0;
            s.vt = sqrt(G*stellar_mass/s.r);

            if (rand()%2 == 1) {
                s.vt *= -1;
            }

            if ( r > 0 ) {
                s.n[0] = (t * num_points) + r - 1;
                s.n[1] = (((t+1)%num_points) * num_points) + r - 1;
                s.n[7] = (((t+num_points-1)%num_points) * num_points) + r - 1;
            }
            if ( r < num_points - 1 ) {
                s.n[3] = (((t+1)%num_points) * num_points) + r + 1;
                s.n[4] = (t * num_points) + r + 1;
                s.n[5] = (((t+num_points-1)%num_points) * num_points) + r + 1;

            }
            s.n[2] = (((t+1)%num_points) * num_points) + r;
            s.n[6] = (((t+num_points-1)%num_points) * num_points) + r;


            segment.push_back(s);
            newSegment.push_back(s);

            SegmentVertices sv;
            sv.v1.x = (radius - radius_step_2) * cos(theta - theta_step_2);
            sv.v1.y = (radius - radius_step_2) * sin(theta - theta_step_2);
            sv.v1.z = 0.0f;
            sv.v1.w = 1.0f;
            sv.v2.x = (radius + radius_step_2) * cos(theta - theta_step_2);
            sv.v2.y = (radius + radius_step_2) * sin(theta - theta_step_2);
            sv.v2.z = 0.0f;
            sv.v2.w = 1.0f;
            sv.v3.x = (radius + radius_step_2) * cos(theta + theta_step_2);
            sv.v3.y = (radius + radius_step_2) * sin(theta + theta_step_2);
            sv.v3.z = 0.0f;
            sv.v3.w = 1.0f;
            sv.v4.x = (radius - radius_step_2) * cos(theta - theta_step_2);
            sv.v4.y = (radius - radius_step_2) * sin(theta - theta_step_2);
            sv.v4.z = 0.0f;
            sv.v4.w = 1.0f;
            sv.v5.x = (radius - radius_step_2) * cos(theta + theta_step_2);
            sv.v5.y = (radius - radius_step_2) * sin(theta + theta_step_2);
            sv.v5.z = 0.0f;
            sv.v5.w = 1.0f;
            sv.v6.x = (radius + radius_step_2) * cos(theta + theta_step_2);
            sv.v6.y = (radius + radius_step_2) * sin(theta + theta_step_2);
            sv.v6.z = 0.0f;
            sv.v6.w = 1.0f;
            segmentVertices.push_back(sv);

            SegmentColours sc;
            sc.c1.r = s.m;
            sc.c1.g = sc.c1.r;
            sc.c1.b = sc.c1.r;
            sc.c1.a = 1.0f;
            sc.c2 = sc.c1;
            sc.c3 = sc.c1;
            sc.c4 = sc.c1;
            sc.c5 = sc.c1;
            sc.c6 = sc.c1;
            segmentColours.push_back(sc);
        }
    }

    fprintf(
            stderr,
            "INFO: segmentVertices.size(): %d\n",
            (unsigned) segmentVertices.size()
    );

    PrintPoints();
    CalcSystemMass();

    InitWindow(argc, argv);

    glewExperimental = GL_TRUE;
    GLenum GlewInitResult = glewInit();

    if (GLEW_OK != GlewInitResult) {
        fprintf(
                stderr,
                "ERROR: %s\n",
                glewGetErrorString(GlewInitResult)
        );
        exit(EXIT_FAILURE);
    }

    fprintf(
            stdout,
            "INFO: OpenGL Version: %s\n",
            glGetString(GL_VERSION)
    );

    CreateShaders();
    CreateVBO();

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

void InitWindow(int argc, char *argv[]) {
    glutInit(&argc, argv);

    glutInitContextVersion(4, 0);
    glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
    glutInitContextProfile(GLUT_CORE_PROFILE);

    glutSetOption(
            GLUT_ACTION_ON_WINDOW_CLOSE,
            GLUT_ACTION_GLUTMAINLOOP_RETURNS
    );

    glutInitWindowSize(CurrentWidth, CurrentHeight);

    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA| GLUT_MULTISAMPLE);

    WindowHandle = glutCreateWindow(WINDOW_TITLE_PREFIX);

    if (WindowHandle < 1) {
        fprintf(
                stderr,
                "ERROR: Could not create a new rendering window.\n"
        );
        exit(EXIT_FAILURE);
    }

    glutReshapeFunc(ResizeFunction);
    glutDisplayFunc(RenderFunction);
    glutIdleFunc(IdleFunction);
    glutTimerFunc(0, TimerFunction, 0);
    glutCloseFunc(Cleanup);
}

void ResizeFunction(int Width, int Height) {
    CurrentWidth = Width;
    CurrentHeight = Height;
    glViewport(0, 0, min(CurrentWidth, CurrentHeight), min(CurrentWidth, CurrentHeight));

}

GLuint programID;
GLuint MatrixID;
glm::mat4 MVP;

void CalcSystemMass() {
    double disk_mass = 0;
    double pr = 0;
    double pt = 0;
    for (size_t i = 0; i < segment.size(); i ++ ) {
        disk_mass += segment[i].m;
        pr += segment[i].m * segment[i].vr;
        pt += segment[i].m * segment[i].vt;
    }
    fprintf(stdout, "INFO: Stellar Mass: %.8e Disk Mass: %.8e : Escape Mass: %.8e Total mass: %.8e, pr = %.8e, pt = %.8e\n",
            stellar_mass, disk_mass, escape_mass, stellar_mass+disk_mass, pr, pt );
}

#define AFAC 0.25f
#define BFAC 0.5f

void ApplyBruteForceGravity() {

    // Calculate acceleration for all points
    for (size_t r = 0; r < num_points; r++) {
        for (size_t t = 0; t < num_points; t++) {

            size_t i1 = (t * num_points) + r;

#ifdef STAR
            newSegment[i1].ar = -1.0 * stellar_mass * G / (segment[i1].r * segment[i1].r);
#else
            newSegment[i1].ar = 0;
#endif
            newSegment[i1].at = 0;

#ifdef SELF
                for (size_t sr = 0; sr < num_points; sr++) {
                    for (size_t st = 0; st < num_points; st++) {

                    size_t i2 = (st * num_points) + sr;

                    if (i2 != i1) {

                        double x = segment[i2].x - segment[i1].x;
                        double y = segment[i2].y - segment[i1].y;
                        double d = sqrt((x * x) * (y * y));
                        double alpha = atan(y/x);
                        double f = d * d;
                        double a = 0;
                        if (f > 0) {
                            a = G * segment[i2].m / f;
                        }
                        double dt = segment[i1].theta - alpha;
                        newSegment[i1].at += a * sin(dt) ; // / (2.0 * M_PI * newSegment[i1].r);
                        newSegment[i1].ar += a * cos(dt) ;
                        //fprintf(1e-22 *  stdout, "INFO: %f, %f, %f, %f, %f, %f\n", x, y, d, alpha, f, a );

                    }
                }
            }
#endif
        }
    }

    // apply acceleration to all velocities
    for ( int i = 0; i < newSegment.size(); i++) {
        newSegment[i].vr += newSegment[i].ar;
        newSegment[i].vt += newSegment[i].at;
        newSegment[i].pr.clear();
        newSegment[i].pt.clear();
    }

    // populate points
    double radius_step = ( outer_radius - inner_radius ) / num_points;
    double radius_step_2 = radius_step / 2;

    // based on velocities, move the matter, conserve momentum
    for (size_t r = 0; r < num_points; r++) {
        for (size_t t = 0; t < num_points; t++) {

            size_t i1 = (t * num_points) + r;

            int * n = &newSegment[i1].n[0];

            double theta_step = 2 * (double)M_PI * newSegment[i1].r / num_points;
            double theta_step_2 = theta_step / 2;
            double r1 = newSegment[i1].r - radius_step_2;
            double r2 = r1 + newSegment[i1].vr;
            double t1 = newSegment[i1].theta - theta_step_2;
            double t2 = t1 + newSegment[i1].vt;

            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            int n1, n2, n3;
            if ( r2 >= r1 ) {
                if ( t2 >= t1 ) {
                    // Q1
                    a1 = min(radius_step, r2-r1)     *min(theta_step, t1+theta_step-t2);
                    a2 = min(radius_step, r2-r1)     *min(theta_step, t2-t1);
                    a3 = min(radius_step, r1+radius_step-r2) *min(theta_step, t2-t1);
                    n1 = 4;
                    n2 = 3;
                    n3 = 2;
                } else {
                    // Q4
                    a1 = min(radius_step, r2-r1)     *min(theta_step, t2+theta_step-t1);
                    a2 = min(radius_step, r2-r1)     *min(theta_step, t1-t2);
                    a3 = min(radius_step, r1+radius_step-r2) *min(theta_step, t1-t2);
                    n1 = 4;
                    n2 = 5;
                    n3 = 6;
                }
            } else {
                if ( t2 >= t1 ) {
                    // Q2
                    a1 = min(radius_step, r1-r2)     *min(theta_step, t1+theta_step-t2);
                    a2 = min(radius_step, r1-r2)     *min(theta_step, t2-t1);
                    a3 = min(radius_step, r2+radius_step-r1) *min(theta_step, t2-t1);
                    n1 = 0;
                    n2 = 1;
                    n3 = 2;
                } else {
                    // Q3
                    a1 = min(radius_step, r1-r2)     *min(theta_step, t2+theta_step-t1);
                    a2 = min(radius_step, r1-r2)     *min(theta_step, t1-t2);
                    a3 = min(radius_step, r2+radius_step-r1) *min(theta_step, t1-t2);
                    n1 = 0;
                    n2 = 7;
                    n3 = 6;
                }

            }

            // area 1
            int i2 = n[n1];
            double dm = newSegment[i1].m*a1/newSegment[i1].area;
            double m1 = newSegment[i1].m;
            if ( i2 < 0 ) {
                if (n1 == 0 ) {
                    stellar_mass += dm;
                } else {
                    escape_mass += dm;
                }
            } else {
                newSegment[i2].pr.push_back({newSegment[i1].vr, dm});
            }
            m1 = max(0.0, m1-dm);

            // area 2
            i2 = n[n2];
            dm = newSegment[i1].m*a2/newSegment[i1].area;
            if ( i2 < 0 ) {
                if (n2 == 1 || n2 == 7 ) {
                    stellar_mass += dm;
                } else {
                    escape_mass += dm;
                }
            } else {

                newSegment[i2].pr.push_back({newSegment[i1].vr, dm/2});
                newSegment[i2].pt.push_back({newSegment[i1].vt, dm/2});

            }
            m1 = max(0.0, m1-dm);


            // area 3
            i2 = n[n3];
            dm = newSegment[i1].m*a3/newSegment[i1].area;
            newSegment[i2].pt.push_back({newSegment[i1].vt, dm});
            m1 = max(0.0, m1-dm);

            newSegment[i1].pr.push_back({newSegment[i1].vr, m1/2});
            newSegment[i1].pt.push_back({newSegment[i1].vt, m1/2});

        }
    }
    for ( int i = 0; i < newSegment.size(); i++) {
        double p = 0.0, m = 0.0;
        for( auto &it : newSegment[i].pt ) {
            p += it.m * it.v;
            m += it.m;
        }
        newSegment[i].vt = 0;
        if ( m > 0 ) {
            newSegment[i].vt = p / m;
        }
        newSegment[i].m = m;
        p = 0.0, m = 0.0;
        for( auto &it : newSegment[i].pr ) {
            p += it.m * it.v;
            m += it.m;
        }
        newSegment[i].vr = 0;
        if ( m > 0 ) {
            newSegment[i].vr = p / m;
        }
        newSegment[i].m += m;
    }

}

void ApplyXNearestNeighbourGravity() {

    std::vector<double> dd(num_points*num_points, 0.0f);
    int span = 3;

    for (size_t t = 0; t < num_points; t++) {

        for (size_t r = 0; r < num_points; r++) {

            size_t i1 = (t * num_points) + r;

            double dr = -1/(segment[i1].r*segment[i1].r);
            double dt = 0;

            for ( int st = -span; st < span; st++) {
                for ( int sr = -span; sr < span; sr++) {

                    size_t tr = max( 0, min((int)num_points, (int)r+sr));
                    size_t tt = ((t+num_points+st)%num_points);
                    size_t i2 = (tt * num_points) + tr;

                    if ( i2 != i1 ) {

                        double x = segment[i2].x - segment[i1].x;
                        double y = segment[i2].y - segment[i1].y;
                        double d = sqrt((x*x)*(y*y));
                        double f = d * d;
                        double nr = 0;
                        if (f > 0) {
                            nr = segment[i2].m / f;
                        }
                        double nt = 0;
                        if ( d > 0) {
                            nt = asin(y / d);
                        }
                        dt += nr * sin(nt);
                        dr += nr * cos(nt);
                    }
                }
            }

            int n[8] = { -1, -1, -1, -1, -1, -1, -1, -1};
            if ( r > 0 ) {
                n[0] = (t * num_points) + r - 1;
                n[1] = (((t+1)%num_points) * num_points) + r - 1;
                n[7] = (((t+num_points-1)%num_points) * num_points) + r - 1;
            }
            if ( r < num_points - 1 ) {
                n[3] = (((t+1)%num_points) * num_points) + r + 1;
                n[4] = (t * num_points) + r + 1;
                n[5] = (((t+num_points-1)%num_points) * num_points) + r + 1;

            }
            n[2] = (((t+1)%num_points) * num_points) + r;
            n[6] = (((t+num_points-1)%num_points) * num_points) + r;


            double d = sqrt(dr*dr+dt*dt);
            d = min( segment[i1].m, d );
            dd[i1] -= d;

            double pdr;
            double pdt;
            double pc = (abs(dr)+abs(dt));
            if ( pc > 0 ) {
                pdr = d * abs(dr) / (abs(dr) + abs(dt));
                pdt = d * abs(dt) / (abs(dr) + abs(dt));
            } else {
                pdr = 0;
                pdt = 0;
            }

            if ( dr > 0 ) {
                if( n[3] >= 0 ) {
                    dd[n[3]] += AFAC * pdr;
                    dd[n[4]] += BFAC * pdr;
                    dd[n[5]] += AFAC * pdr;
                }
            } else {
                if( n[1] >= 0 ) {
                    dd[n[1]] += AFAC * pdr;
                    dd[n[0]] += BFAC * pdr;
                    dd[n[7]] += AFAC * pdr;
                }
            }

            if ( dt > 0 ) {
                dd[n[2]] += BFAC * pdt;
                if( n[1] >= 0 ) {

                    dd[n[1]] += AFAC * pdt;
                }
                if( n[3] >= 0 ) {
                    dd[n[3]] += AFAC * pdt;
                }
            } else {
                if( n[7] >= 0 ) {
                    dd[n[7]] += AFAC * pdt;
                }
                    if( n[5] >= 0 ) {
                        dd[n[5]] += AFAC * pdt;
                    }
                dd[n[6]] += BFAC * pdt;
            }

        }
    }

    for ( int i = 0; i < newSegment.size(); i++) {
        newSegment[i].m = max( 0.0, newSegment[i].m + 1e16*dd[i]);
    }
}



double r2 = 1/sqrt(2);


void ApplyNearestNeighbourGravity() {

    std::vector<double> dd(num_points*num_points, 0.0);


    for (size_t t = 0; t < num_points; t++) {

        for (size_t r = 1; r < num_points-1; r++) {

            size_t i = (t * num_points) + r;
            size_t n[8];
            n[0] = (t * num_points) + r - 1;
            n[1] = (((t+1)%num_points) * num_points) + r - 1;
            n[2] = (((t+1)%num_points) * num_points) + r;
            n[3] = (((t+1)%num_points) * num_points) + r + 1;
            n[4] = (t * num_points) + r + 1;
            n[5] = (((t+num_points-1)%num_points) * num_points) + r + 1;
            n[6] = (((t+num_points-1)%num_points) * num_points) + r;
            n[7] = (((t+num_points-1)%num_points) * num_points) + r - 1;

            // radial component
            double dr =   ( r2 * segment[n[3]].m + segment[n[4]].m + r2 * segment[n[5]].m )
                       - ( r2 * segment[n[1]].m + segment[n[0]].m + r2 * segment[n[7]].m )
                       - 10/(segment[i].r*segment[i].r);

            // angular component
            double dt =   ( r2 * segment[n[1]].m + segment[n[2]].m + r2 * segment[n[3]].m )
                       - ( r2 * segment[n[5]].m + segment[n[6]].m + r2 * segment[n[7]].m );

            double d = sqrt(dr*dr+dt*dt);
            d = min( segment[i].m, d );
            dd[i] -= d;

            double pdr;
            double pdt;
            double pc = (abs(dr)+abs(dt));
            if ( pc > 0 ) {
                pdr = d * abs(dr) / (abs(dr) + abs(dt));
                pdt = d * abs(dt) / (abs(dr) + abs(dt));
            } else {
                pdr = 0;
                pdt = 0;
            }
            if ( dr > 0 ) {
                dd[n[3]] += 0.25 * pdr;
                dd[n[4]] += 0.50 * pdr;
                dd[n[5]] += 0.25 * pdr;
            } else {
                dd[n[1]] += 0.25 * pdr;
                dd[n[0]] += 0.50 * pdr;
                dd[n[7]] += 0.25 * pdr;
            }

            if ( dt > 0 ) {
                dd[n[1]] += 0.25 * pdt;
                dd[n[2]] += 0.50 * pdt;
                dd[n[3]] += 0.25 * pdt;
            } else {
                dd[n[7]] += 0.25 * pdt;
                dd[n[6]] += 0.50 * pdt;
                dd[n[5]] += 0.25 * pdt;
            }

        }
    }

    for ( int i = 0; i < newSegment.size(); i++) {
        newSegment[i].m = max( 0.0, newSegment[i].m + dd[i]);
    }
}

void RenderFunction(void) {
    ++FrameCount;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Use our shader
    glUseProgram(programID);

    // Send our transformation to the currently bound shader,
    // in the "MVP" uniform
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);


    glDrawArrays(GL_TRIANGLES, 0, segmentVertices.size()*6);

    glFlush();
    glutSwapBuffers();

    ApplyBruteForceGravity();

    segment = newSegment;

    CalcSystemMass();

    MapSegmentToColor();

    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentVertices.size() * sizeof(SegmentVertices), &segmentVertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(XY), 0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentColours.size() * sizeof(SegmentColours), &segmentColours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(RGBA), 0);
    glEnableVertexAttribArray(0);

}

void MapSegmentToColor() {

    double d_max = 250;

    for (int i = 0; i < segmentColours.size(); i++) {
        SegmentColours* scp = &segmentColours[i];
        scp->c1.r = min((float)((segment[i].m/segment[i].area)/250), 1.0f);
        scp->c1.g = scp->c1.r;
        scp->c1.b = scp->c1.r;
        scp->c2 = scp->c1;
        scp->c3 = scp->c1;
        scp->c4 = scp->c1;
        scp->c5 = scp->c1;
        scp->c6 = scp->c1;
    }
}

void IdleFunction(void) {
    glutPostRedisplay();
}

void TimerFunction(int Value) {
    if (0 != Value) {
        char *TempString = (char *)
                malloc(512 + strlen(WINDOW_TITLE_PREFIX));

        sprintf(
                TempString,
                "%s: %d Frames Per Second @ %d x %d",
                WINDOW_TITLE_PREFIX,
                FrameCount * 4,
                CurrentWidth,
                CurrentHeight
        );

        glutSetWindowTitle(TempString);
        free(TempString);
    }

    FrameCount = 0;
    glutTimerFunc(250, TimerFunction, 1);
}

void Cleanup(void) {
    DestroyShaders();
    DestroyVBO();
}

void CreateVBO(void) {
    GLenum ErrorCheckValue = glGetError();

    glGenVertexArrays(1, &VaoId);
    glBindVertexArray(VaoId);

    glGenBuffers(1, &VboAreaId);
    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentVertices.size() * sizeof(SegmentVertices), &segmentVertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(XY), 0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &VboPropertyId);
    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentColours.size() * sizeof(SegmentColours), &segmentColours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(RGBA), 0);
    glEnableVertexAttribArray(1);

    ErrorCheckValue = glGetError();
    if (ErrorCheckValue != GL_NO_ERROR) {
        fprintf(
                stderr,
                "ERROR: Could not create a VBO: %s \n",
                gluErrorString(ErrorCheckValue)
        );

        exit(-1);
    }
}

void DestroyVBO(void) {
    GLenum ErrorCheckValue = glGetError();

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDeleteBuffers(1, &VboPropertyId);

    glBindVertexArray(0);
    glDeleteVertexArrays(1, &VaoId);

    ErrorCheckValue = glGetError();
    if (ErrorCheckValue != GL_NO_ERROR) {
        fprintf(
                stderr,
                "ERROR: Could not destroy the VBO: %s \n",
                gluErrorString(ErrorCheckValue)
        );

        exit(-1);
    }
}

void CreateShaders(void) {
    // Create and compile our GLSL program from the shaders
    programID = LoadShaders("../SimpleTransform.vertexshader",
                            "../SingleColor.fragmentshader");

    // Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    //glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
    // Or, for an ortho camera :
    float ortho = (float)outer_radius;
    float factor = 1.1f;

    glm::mat4 Projection = glm::ortho(-1.0f * factor * ortho, factor * ortho, -1.0f * factor * ortho, factor * ortho, -1.0f, 1.0f); // In world coordinates
    fprintf(
            stderr,
            "INFO: ortho: %f\n",
            ortho
    );


    // Camera matrix
    glm::mat4 View = glm::lookAt(
            glm::vec3(0, 0, 0.1), // Camera is at (4,3,3), in World Space
            glm::vec3(0, 0, 0), // and looks at the origin
            glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
    );
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    MVP = Projection * View * Model; // Remember, matrix multiplication is the other way around

    // Get a handle for our "MVP" uniform
    MatrixID = glGetUniformLocation(programID, "MVP");


}

void DestroyShaders(void) {
    glDeleteProgram(programID);
}

void PrintPoints() {
    for (int i = 0; i < segmentVertices.size(); i++) {
        SegmentVertices pt = segmentVertices[i];
        fprintf(stdout,
                "%d:([%f,%f],[%f,%f],[%f,%f]) ",
                i, pt.v1.x, pt.v1.y, pt.v2.x, pt.v2.y, pt.v3.x, pt.v3.y);
    }
    fprintf(stdout, "\n");
}