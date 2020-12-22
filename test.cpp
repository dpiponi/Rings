#include <math.h>
#include <stdio.h>
#include <thread>
#include <vector>

#include "raylib.h"

#define n_particles 100000
#define n_threads 8

double M = 1.0;
double G = 1.0;
double Distance = 1.0;
double AngularVelocity;
double MoonMass = 0.000001;
double MoonRadius = 1.125;

double dt = 0.0001;

Color ParticleColour = {255, 255, 255, 60};
Color MoonColour = {255, 255, 255, 255};

struct V2d
{
  double x;
  double y;

  V2d& operator+=(const V2d& Other)
  {
    x += Other.x;
    y += Other.y;

    return *this;
  }

  V2d operator-(const V2d& Other)
  {
    return V2d{x - Other.x, y - Other.y};
  }

  V2d operator+(const V2d& Other)
  {
    return V2d{x + Other.x, y + Other.y};
  }
};

V2d operator*(double A, const V2d& V)
{
  return V2d{A * V.x, A * V.y};
}

V2d particles[n_particles];
V2d velocities[n_particles];

V2d Shepherd;

struct M22d
{
  double M[2][2];
};

M22d Rotation(double angle)
{
  double C = cos(angle);
  double S = sin(angle);
  return M22d{C, -S, S, C};
};

double OrbitalVelocityFromRadius(double Radius)
{
  return sqrt(G * M / Radius);
}

void InitDynamics()
{
  AngularVelocity = sqrt(G * M / (Distance * Distance * Distance));

  for (int i = 0; i < n_particles; ++i)
  {
    double Radius = 1.0 + 0.1 * 0.0001 * GetRandomValue(0, 10000);
    double Angle = 2 * M_PI * 0.000001 * GetRandomValue(0, 1000000);
    particles[i].x = Radius * cos(Angle);
    particles[i].y = Radius * sin(Angle);
    double Speed = OrbitalVelocityFromRadius(Radius);
    velocities[i].x = particles[i].y / Radius * Speed;
    velocities[i].y = -particles[i].x / Radius * Speed;
  }
}

// Returns acceleration
V2d Gravitation(double M, const V2d& Position)
{
  double R = hypot(Position.x, Position.y);
  double ScalarF = G * M / (R * R);
  double FX = -ScalarF * Position.x / R;
  double FY = -ScalarF * Position.y / R;
  return V2d{FX, FY};
}

V2d Shepherd0, Shepherd1;

double t = 0.;
void DynamicsRange(int From, int To, int NSteps)
{
  for (int Step = 0; Step < NSteps; ++Step)
  {
    // XXX Combine
    for (int i = From; i < To; ++i)
    {
      velocities[i] += 0.5 * dt * (Gravitation(M, particles[i]) + Gravitation(MoonMass, particles[i] - Shepherd0));
      particles[i] += dt * velocities[i];
      velocities[i] += 0.5 * dt * (Gravitation(M, particles[i]) + Gravitation(MoonMass, particles[i] - Shepherd1));
    }

    t += dt;
  }
  Shepherd = Shepherd1; // XXX
}

void dynamics()
{
    double AngularVelocity = OrbitalVelocityFromRadius(MoonRadius) / MoonRadius;
    double angle0 = AngularVelocity * t;
    double angle1 = AngularVelocity * (t + 0.5 * dt);
//     printf("AngularVelocity=%f\n", AngularVelocity);
//     printf("angle0=%f\n", angle0);
    Shepherd0 = MoonRadius * V2d{cos(angle0), -sin(angle0)};
    Shepherd1 = MoonRadius * V2d{cos(angle1), -sin(angle1)};

  std::vector<std::thread> Threads;
  for (int Thread = 0; Thread < n_threads; ++Thread)
  {
    Threads.emplace_back([Thread]()
    {
      int First = Thread * n_particles / n_threads;
      int Last = std::min((Thread + 1) * n_particles / n_threads, n_particles);
      DynamicsRange(First, Last, 100);
    });
  }

  for (int Thread = 0; Thread < n_threads; ++Thread)
  {
    Threads[Thread].join();
  }

//   printf("%f %f\n", particles[0].x, particles[0].y);
}

const int ScreenWidth = 1400;
const int ScreenHeight = 1400;

void DrawParticles()
{
#if 0
  double Top = 1.2;
  double Bottom = 0.8;
  double Left = -0.1;
  double Right = 0.1;
#else
  double Top = 1.2;
  double Bottom = -1.2;
  double Left = -1.2;
  double Right = 1.2;
#endif

  double Distance = MoonRadius;
  double AngularVelocity = sqrt(G * M / (Distance * Distance * Distance));
  M22d Rotate = Rotation(-t * AngularVelocity);

  for (int i = 0; i < n_particles; ++i)
  {
    V2d P = particles[i];
    if (P.x > Left && P.x < Right &&
        P.y > Bottom && P.y < Top)
    {
      int ScreenX = rint((P.x - Left) / (Right - Left) * ScreenWidth);
      int ScreenY = rint((P.y - Top) / (Bottom - Top) * ScreenHeight);
      DrawPixel(ScreenX, ScreenY, ParticleColour);
    }
  }

  int ScreenX = rint((Shepherd.x - Left) / (Right - Left) * ScreenWidth);
  int ScreenY = rint((Shepherd.y - Top) / (Bottom - Top) * ScreenHeight);
  DrawCircle(ScreenX, ScreenY, 4.0, MoonColour);
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
