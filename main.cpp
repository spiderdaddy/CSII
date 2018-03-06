
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

struct Segment {
    float r, theta, x, y;
    float area;
    float density;
};


int num_points = 900;
float inner_radius = .39;
float outer_radius = 40;

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

void MapDensityToColour();

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
    float theta_step = 2 * M_PI / num_points;
    float theta_step_2 = theta_step / 2;
    float radius_step = ( outer_radius - inner_radius ) / num_points;
    float radius_step_2 = radius_step / 2;

    for (size_t t = 0; t < num_points; t++) {

        float theta = t * theta_step;

        for (size_t r = 0; r < num_points; r++) {

            float radius = inner_radius + (r * radius_step) + radius_step_2;

            Segment s;
            s.r = radius;
            s.theta = theta;
            s.x = radius * cos(theta);
            s.y = radius * sin(theta);
            s.area = 1;
            s.density = float(rand() % 64) / 255;
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
            sc.c1.r = s.density;
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

void CalcTotalMass() {
    float mass = 0;
    for (size_t i = 0; i < segment.size(); i ++ ) {
        mass += segment[i].density * segment[i].area;
    }
    fprintf(stdout, "INFO: Total mass: %f\n", mass );
}

void ApplyStellarGravity() {


    for (size_t t = 0; t < num_points; t++) {

        for (size_t r = 1; r < num_points-1; r++) {

            size_t i = (t * num_points) + r;
            size_t n[8];
            n[0] = (t * num_points) + r - 1;
            n[1] = (((t+1)%num_points) * num_points) + r - 1;
            n[2] = (((t+1)%num_points) * num_points) + r;
            n[3] = (((t+1)%num_points) * num_points) + r + 1;
            n[4] = (t * num_points) + r + 1;
            n[5] = (((t+num_points-1)%num_points) * num_points) + r - 1;
            n[6] = (((t+num_points-1)%num_points) * num_points) + r;
            n[7] = (((t+num_points-1)%num_points) * num_points) + r + 1;

            // radial component
            float d = ( 0.707 * segment[n[3]].density + segment[n[4]].density + 0.707 * segment[n[5]].density )
                    - ( 0.707 * segment[n[1]].density + segment[n[0]].density + 0.707 * segment[n[7]].density );

            newSegment[i].density -= d;
            if ( d > 0 ) {
                newSegment[n[4]].density += d;
            } else {
                newSegment[n[0]].density += d;
            }

            // angular component
            d = ( 0.707 * segment[n[3]].density + segment[n[2]].density + 0.707 * segment[n[5]].density )
                      - ( 0.707 * segment[n[1]].density + segment[n[6]].density + 0.707 * segment[n[7]].density );

            newSegment[i].density -= d;
            if ( d > 0 ) {
                newSegment[n[6]].density += d;
            } else {
                newSegment[n[2]].density += d;
            }
        }
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

    ApplyStellarGravity();

    segment = newSegment;

    MapDensityToColour();

    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentVertices.size() * sizeof(SegmentVertices), &segmentVertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(XY), 0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentColours.size() * sizeof(SegmentColours), &segmentColours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(RGBA), 0);
    glEnableVertexAttribArray(0);

    CalcTotalMass();

}

void MapDensityToColour() {
    for (int i = 0; i < segmentColours.size(); i++) {
        SegmentColours* scp = &segmentColours[i];
        scp->c1.r = segment[i].density;
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
    float ortho = outer_radius;
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