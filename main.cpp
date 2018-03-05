
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

struct SegmentArea {
    XY coord1;
    XY coord2;
    XY coord3;
    XY coord4;
    XY coord5;
    XY coord6;

};

struct SegmentProperty {
    RGBA colour1;
    RGBA colour2;
    RGBA colour3;
    RGBA colour4;
    RGBA colour5;
    RGBA colour6;

};

std::vector<SegmentArea> segmentAreas;
std::vector<SegmentProperty> segmentProps;


int num_points = 300;
float inner_radius = .39;
float outer_radius = 40;

int
        CurrentWidth = 800,
        CurrentHeight = 800,
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

int main(int argc, char *argv[]) {
    Initialize(argc, argv);

    glutMainLoop();

    exit(EXIT_SUCCESS);
}


void Initialize(int argc, char *argv[]) {

    // populate points
    float theta_step = 2 * M_PI / num_points;
    float theta_step_2 = theta_step / 2;
    float radius_step = ( outer_radius - inner_radius ) / num_points;
    float radius_step_2 = radius_step / 2;
    for (size_t r = 0; r < num_points; r++) {

        float radius = inner_radius + (r * radius_step) + radius_step_2;

        for (size_t t = 0; t < num_points; t++) {

            float theta = t * theta_step;

            SegmentArea pt;
            pt.coord1.x = (radius - radius_step_2) * cos(theta - theta_step_2);
            pt.coord1.y = (radius - radius_step_2) * sin(theta - theta_step_2);
            pt.coord1.z = 0.0f;
            pt.coord1.w = 1.0f;
            pt.coord2.x = (radius + radius_step_2) * cos(theta - theta_step_2);
            pt.coord2.y = (radius + radius_step_2) * sin(theta - theta_step_2);
            pt.coord2.z = 0.0f;
            pt.coord2.w = 1.0f;
            pt.coord3.x = (radius + radius_step_2) * cos(theta + theta_step_2);
            pt.coord3.y = (radius + radius_step_2) * sin(theta + theta_step_2);
            pt.coord3.z = 0.0f;
            pt.coord3.w = 1.0f;
            pt.coord4.x = (radius - radius_step_2) * cos(theta - theta_step_2);
            pt.coord4.y = (radius - radius_step_2) * sin(theta - theta_step_2);
            pt.coord4.z = 0.0f;
            pt.coord4.w = 1.0f;
            pt.coord5.x = (radius - radius_step_2) * cos(theta + theta_step_2);
            pt.coord5.y = (radius - radius_step_2) * sin(theta + theta_step_2);
            pt.coord5.z = 0.0f;
            pt.coord5.w = 1.0f;
            pt.coord6.x = (radius + radius_step_2) * cos(theta + theta_step_2);
            pt.coord6.y = (radius + radius_step_2) * sin(theta + theta_step_2);
            pt.coord6.z = 0.0f;
            pt.coord6.w = 1.0f;

            segmentAreas.push_back(pt);

            SegmentProperty sp;
            sp.colour1.r = float(rand() % 256) / 255;
            sp.colour1.g = sp.colour1.r;
            sp.colour1.b = sp.colour1.r;
            sp.colour1.a = 1.0f;
            sp.colour2 = sp.colour1;
            sp.colour3 = sp.colour1;
            sp.colour4 = sp.colour1;
            sp.colour5 = sp.colour1;
            sp.colour6 = sp.colour1;
            segmentProps.push_back(sp);
        }
    }

    fprintf(
            stderr,
            "INFO: segmentAreas.size(): %d\n",
            (unsigned) segmentAreas.size()
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


void RenderFunction(void) {
    ++FrameCount;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Use our shader
    glUseProgram(programID);

    // Send our transformation to the currently bound shader,
    // in the "MVP" uniform
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);


    glDrawArrays(GL_TRIANGLES, 0, segmentAreas.size()*6);

    glFlush();
    glutSwapBuffers();

    SegmentProperty pt = segmentProps[segmentProps.size() - 1];
    pt.colour1.r = segmentProps[0].colour1.r;
    pt.colour1.g = segmentProps[0].colour1.g;
    pt.colour1.b = segmentProps[0].colour1.b;
    pt.colour2 = segmentProps[0].colour1;
    pt.colour3 = segmentProps[0].colour1;
    pt.colour4 = segmentProps[0].colour1;
    pt.colour5 = segmentProps[0].colour1;
    pt.colour6 = segmentProps[0].colour1;
    segmentProps[segmentProps.size() - 1] = pt;

    for (int i = 0; i < segmentProps.size() - 1; i++) {
        SegmentProperty pt = segmentProps[i];
        pt.colour1.r = segmentProps[i + 1].colour1.r;
        pt.colour1.g = segmentProps[i + 1].colour1.g;
        pt.colour1.b = segmentProps[i + 1].colour1.b;
        pt.colour2 = segmentProps[i + 1].colour1;
        pt.colour3 = segmentProps[i + 1].colour1;
        pt.colour4 = segmentProps[i + 1].colour1;
        pt.colour5 = segmentProps[i + 1].colour1;
        pt.colour6 = segmentProps[i + 1].colour1;
        segmentProps[i] = pt;
    }

    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentAreas.size() * sizeof(SegmentArea), &segmentAreas[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(XY), 0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentProps.size() * sizeof(SegmentProperty), &segmentProps[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(RGBA), 0);
    glEnableVertexAttribArray(0);



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
    glBufferData(GL_ARRAY_BUFFER, segmentAreas.size() * sizeof(SegmentArea), &segmentAreas[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(XY), 0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &VboPropertyId);
    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentProps.size() * sizeof(SegmentProperty), &segmentProps[0], GL_STATIC_DRAW);
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
    for (int i = 0; i < segmentAreas.size(); i++) {
        SegmentArea pt = segmentAreas[i];
        fprintf(stdout,
                "%d:([%f,%f],[%f,%f],[%f,%f]) ",
                i, pt.coord1.x, pt.coord1.y, pt.coord2.x, pt.coord2.y, pt.coord3.x, pt.coord3.y);
    }
    fprintf(stdout, "\n");
}