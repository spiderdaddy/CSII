//
// Created by lf on 20/03/18.
//

#include <cstdio>
#include <vector>
#include <cstring>

#include <GL/glew.h>
#include <GL/freeglut.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <shader.hpp>

#include "graphics.h"

#define WINDOW_TITLE_PREFIX "Gravity"

using namespace glm;

int CurrentWidth = 1800;
int CurrentHeight = 900;
int WindowHandle = 0;
unsigned FrameCount = 0;

GLuint VaoId;
GLuint VboAreaId;
GLuint VboPropertyId;

GLuint programID;
GLuint MatrixID;

glm::mat4 MVP;

void GraphicsMainLoop() {

    glutMainLoop();

}

Disk * disk;

void InitializeGraphics(int argc, char *argv[], Disk *d) {

    disk = d;

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


void RenderFunction(void) {

    ++FrameCount;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Use our shader
    glUseProgram(programID);

    // Send our transformation to the currently bound shader,
    // in the "MVP" uniform
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

    std::vector<Disk::SegmentVertices> segmentVertices = disk->getSegmentVertices();
    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentVertices.size() * sizeof(Disk::SegmentVertices), &segmentVertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Disk::XYZW_GL), 0);
    glEnableVertexAttribArray(1);

    disk->MapSegmentToColor();
    std::vector<Disk::SegmentColours> segmentColours = disk->getSegmentColours();
    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentColours.size() * sizeof(Disk::SegmentColours), &segmentColours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Disk::RGBA_GL), 0);
    glEnableVertexAttribArray(0);

    glDrawArrays(GL_TRIANGLES, 0, segmentVertices.size()*6);

    std::vector<Disk::SegmentVertices> gravityVertices = disk->getGravityVertices();
    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentVertices.size() * sizeof(Disk::SegmentVertices), &gravityVertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Disk::XYZW_GL), 0);
    glEnableVertexAttribArray(1);

    disk->MapGravityToColor();
    std::vector<Disk::SegmentColours> gravityColours = disk->getGravityColours();
    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentColours.size() * sizeof(Disk::SegmentColours), &gravityColours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Disk::RGBA_GL), 0);
    glEnableVertexAttribArray(0);

    glDrawArrays(GL_TRIANGLES, 0, gravityVertices.size()*6);


    glFlush();
    glutSwapBuffers();

    ApplyGravity(disk, 6, "test");

    disk->swapSegments();

    disk->CalcSystemMass();
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
    unsigned int ErrorCheckValue = glGetError();

    std::vector<Disk::SegmentVertices> segmentVertices = disk->getSegmentVertices();
    std::vector<Disk::SegmentColours> segmentColours = disk->getSegmentColours();

    glGenVertexArrays(1, &VaoId);
    glBindVertexArray(VaoId);

    glGenBuffers(1, &VboAreaId);
    glBindBuffer(GL_ARRAY_BUFFER, VboAreaId);
    glBufferData(GL_ARRAY_BUFFER, segmentVertices.size() * sizeof(Disk::SegmentVertices), &segmentVertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Disk::XYZW_GL), 0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &VboPropertyId);
    glBindBuffer(GL_ARRAY_BUFFER, VboPropertyId);
    glBufferData(GL_ARRAY_BUFFER, segmentColours.size() * sizeof(Disk::SegmentColours), &segmentColours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Disk::RGBA_GL), 0);
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
    unsigned int ErrorCheckValue = glGetError();

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
    // glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
    // Or, for an ortho camera :
    float xortho = (float) OUTER_RADIUS * 2;
    float yortho = (float) OUTER_RADIUS;
    float factor = 1.1f;

    mat4x4 Projection = glm::ortho(-1.0f * factor * xortho, factor * xortho, -1.0f * factor * yortho, factor * yortho, -2.0f, 2.0f); // In world coordinates
    fprintf(
            stderr,
            "INFO: xortho: %f yortho: %f\n",
            xortho, yortho
    );

    // Camera matrix
    mat4x4 View = glm::lookAt(
            highp_vec3(0, 0, 1), // Camera in World Space
            highp_vec3(0, 0, 0), // and looks at the origin
            highp_vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
    );
    // Model matrix : an identity matrix (model will be at the origin)
    mat4x4 Model = mat4x4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    MVP = Projection * View * Model; // Remember, matrix multiplication is the other way around

    // Get a handle for our "MVP" uniform
    MatrixID = glGetUniformLocation(programID, "MVP");


}

void DestroyShaders(void) {
    glDeleteProgram(programID);
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
    //glViewport(0, 0, std::min(CurrentWidth, CurrentHeight), std::min(CurrentWidth, CurrentHeight));
    glViewport(0, 0, CurrentWidth, CurrentHeight);

}

