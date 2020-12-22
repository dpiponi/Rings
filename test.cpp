#include <math.h>
#include <stdio.h>
#include <thread>

#include "raylib.h"

#define n_particles 100000
#define n_threads 8

double M = 1.0;
double G = 1.0;
double Distance = 1.0;
double AngularVelocity;

double dt = 0.0001;

// Color ParticleColour = {255, 255, 255, 25};
Color ParticleColour = {255, 255, 255, 45};

typedef struct
{
  double x;
  double y;
} V2d;

V2d particles[n_particles];
V2d velocities[n_particles];

double OrbitalVelocityFromRadius(double Radius)
{
  return sqrt(G * M * Radius);
}

void InitDynamics()
{
  AngularVelocity = sqrt(G * M / (Distance * Distance * Distance));

  for (int i = 0; i < n_particles; ++i)
  {
    double Radius = 0.9 + 0.2 * 0.0001 * GetRandomValue(0, 10000);
    double Angle = 2 * M_PI * 0.000001 * GetRandomValue(0, 1000000);
    particles[i].x = Radius * cos(Angle);
    particles[i].y = Radius * sin(Angle);
    double Speed = 0.5*OrbitalVelocityFromRadius(Radius);
    velocities[i].x = particles[i].y / Radius * Speed;
    velocities[i].y = -particles[i].x / Radius * Speed;
  }
}

void DynamicsRange(int From, int To, int NSteps)
{
  for (int Step = 0; Step < NSteps; ++Step)
  {
    for (int i = From; i < To; ++i)
    {
      double R = hypot(particles[i].x, particles[i].y);
      double ScalarF = G * M / (R * R);
      double FX = -ScalarF * particles[i].x / R;
      double FY = -ScalarF * particles[i].y / R;
      velocities[i].x += dt * FX;
      velocities[i].y += dt * FY;
      particles[i].x += dt * velocities[i].x;
      particles[i].y += dt * velocities[i].y;
    }
  }
}

void dynamics()
{
  printf("Dynamics\n");
  {
    DynamicsRange(0, n_particles, 10);
  }

  printf("%f %f\n", particles[0].x, particles[0].y);
}

const int ScreenWidth = 1000;
const int ScreenHeight = 1000;

void DrawParticles()
{
#if 0
  double Top = 1.01;
  double Bottom = 0.99;
  double Left = -0.01;
  double Right = 0.01;
#else
  double Top = 1.;
  double Bottom = -1.;
  double Left = -1.;
  double Right = 1.;
#endif

  for (int i = 0; i < n_particles; ++i)
  {
    if (particles[i].x > Left && particles[i].x < Right &&
        particles[i].y > Bottom && particles[i].y < Top)
    {
      int ScreenX = rint((particles[i].x - Left) / (Right - Left) * ScreenWidth);
      int ScreenY = rint((particles[i].y - Top) / (Bottom - Top) * ScreenHeight);
      DrawPixel(ScreenX, ScreenY, ParticleColour);
    }
  }
}

int main(void)
{

    InitDynamics();

    InitWindow(ScreenWidth, ScreenHeight, "raylib [core] example - basic window");

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second
    
    float t = 0.0;
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
      t += 0.001;
      dynamics();

      BeginDrawing();
      Color Black = {0, 0, 0, 255};
      ClearBackground(Black);
      DrawParticles();
      EndDrawing();
    }
    CloseWindow();        // Close window and OpenGL context

    return 0;
}
