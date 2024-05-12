#include "pr.hpp"
#include <iostream>
#include <memory>

template <class T>
inline T ss_max(T x, T y)
{
    return (x < y) ? y : x;
}

template <class T>
inline T ss_min(T x, T y)
{
    return (y < x) ? y : x;
}

float sign_of(float v)
{
    return v < 0.0f ? -1.0f : 1.0f;
}

// ax^2 + bx + c == 0
int solve_quadratic(float xs[2], float a, float b, float c)
{
    float det = b * b - 4.0f * a * c;
    if (det < 0.0f)
    {
        return 0;
    }

    float k = (-b - sign_of(b) * std::sqrtf(det)) / 2.0f;
    float x0 = k / a;
    float x1 = c / k;
    xs[0] = ss_min(x0, x1);
    xs[1] = ss_max(x0, x1);
    return 2;
}
auto sqr = [](float x) { return x * x; };

float intersect_ray_ellipsoid(glm::vec3 U, glm::vec3 V, glm::vec3 W, glm::vec3 ro, glm::vec3 rd)
{
    glm::vec3 u = U / glm::dot(U, U);
    glm::vec3 v = V / glm::dot(V, V);
    glm::vec3 w = W / glm::dot(W, W);

    float k = 1.0f;
    float t_delta = -glm::dot(rd, ro) / glm::dot(rd, rd);
    glm::vec3 ro_prime = ro + rd * t_delta;

    float urd = glm::dot(u, rd);
    float vrd = glm::dot(v, rd);
    float wrd = glm::dot(w, rd);
    float uro = glm::dot(u, ro_prime);
    float vro = glm::dot(v, ro_prime);
    float wro = glm::dot(w, ro_prime);
    float A = sqr(urd) + sqr(vrd) + sqr(wrd);
    float B = 2.0f * (urd * uro + vrd * vro + wrd * wro);
    float C = sqr(uro) + sqr(vro) + sqr(wro) - k * k;

    float xs[2];
    if (solve_quadratic(xs, A, B, C))
    {
        return xs[0] + t_delta;
    }
    return -1.0f;
}

float lengthSquared(glm::vec3 v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
// lambda0 is larger
void eignValues(float* lambda0, float* lambda1, float* determinant, const glm::mat2& mat)
{
    float mean = (mat[0][0] + mat[1][1]) * 0.5f;
    float det = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    float d = std::sqrtf(ss_max(mean * mean - det, 0.0f));
    *lambda0 = mean + d;
    *lambda1 = mean - d;
    *determinant = det;
}

void eigenVectors_of_symmetric(glm::vec2* eigen0, glm::vec2* eigen1, const glm::mat2& m, float lambda)
{
    float s11 = m[0][0];
    float s22 = m[1][1];
    float s12 = m[1][0];

    // to workaround lambda0 == lambda1
    float eps = 1e-15f;
    glm::vec2 e0 = glm::normalize(s11 < s22 ? glm::vec2(s12 + eps, lambda - s11) : glm::vec2(lambda - s22, s12 + eps));
    glm::vec2 e1 = { -e0.y, e0.x };
    *eigen0 = e0;
    *eigen1 = e1;
}

int main() {
    using namespace pr;

    SetDataDir(ExecutableDir());

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 4, 4, 4 };
    camera.lookat = { 0, 0, 0 };
    camera.zNear = 1.0f;
    camera.zFar = 1000.0f;

    //camera.origin *= 100.0f;
    //camera.fovy = 0.005f;

    double e = GetElapsedTime();

    ITexture* tex = CreateTexture();
    Image2DRGBA8 image;
    int stride = 1;

    Image2DRGBA8 earth;
    earth.load("earth.jpg");

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        ClearBackground(0.1f, 0.1f, 0.1f, 1);

        BeginCamera(camera);

        PushGraphicState();

        DrawGrid(GridAxis::XZ, 1.0f, 10, { 128, 128, 128 });

        static glm::vec3 U = { 1, 0, 0 };
        static glm::vec3 V = { 0, 1, 0 };
        static glm::vec3 W = { 0, 0, 1 };
        ManipulatePosition(camera, &U, 0.2f);
        ManipulatePosition(camera, &V, 0.2f);
        ManipulatePosition(camera, &W, 0.2f);

        {
            glm::vec3 N = glm::cross(U, V);
            V = glm::normalize(glm::cross(N, U)) * glm::length(V);

            DrawArrow({ 0,0,0 }, U, 0.01f, { 255,0,0 });
            DrawArrow({ 0,0,0 }, V, 0.01f, { 0,255,0 });
            DrawArrow({ 0,0,0 }, W, 0.01f, { 0,0,255 });

            W = glm::normalize(N) * glm::length(W);
        }


        // https://tavianator.com/2014/ellipsoid_bounding_boxes.html
        float dx = std::sqrt(U.x * U.x + V.x * V.x + W.x * W.x);
        float dy = std::sqrt(U.y * U.y + V.y * V.y + W.y * W.y);
        float dz = std::sqrt(U.z * U.z + V.z * V.z + W.z * W.z);
        DrawCube({ 0, 0, 0 }, { dx * 2, dy * 2, dz * 2 }, { 255,255,255 });

        
        glm::mat3 R = glm::mat3(glm::normalize(U), glm::normalize(V), glm::normalize(W));
        glm::mat3 L = glm::mat3(
            lengthSquared(U), 0, 0,
            0, lengthSquared(V), 0,
            0, 0, lengthSquared(W)
        );
        glm::mat3 cov = R * L * glm::transpose(R);
        glm::mat3 inv_cov = glm::inverse( cov );
        
        { // xy projection

            // need to take 2x2 in cov not inv_cov
            glm::mat2 cov_pj = glm::mat2(cov[0][0], cov[1][0], cov[0][1], cov[1][1]);
            float det;
            float lambda0;
            float lambda1;
            eignValues(&lambda0, &lambda1, &det, cov_pj);

            glm::vec2 e0;
            glm::vec2 e1;
            eigenVectors_of_symmetric(&e0, &e1, cov_pj, lambda0);
            e0 *= std::sqrtf(lambda0);
            e1 *= std::sqrtf(lambda1);

            DrawEllipse(glm::vec3(0, 0, dz), glm::vec3(e0.x, e0.y, 0), glm::vec3(e1.x, e1.y, 0), { 128, 0, 0 });
        }

        { // yz projection

            // need to take 2x2 in cov not inv_cov
            glm::mat2 cov_pj = glm::mat2(cov[1][1], cov[2][1], cov[1][2], cov[2][2]);
            float det;
            float lambda0;
            float lambda1;
            eignValues(&lambda0, &lambda1, &det, cov_pj);

            glm::vec2 e0;
            glm::vec2 e1;
            eigenVectors_of_symmetric(&e0, &e1, cov_pj, lambda0);
            e0 *= std::sqrtf(lambda0);
            e1 *= std::sqrtf(lambda1);

            DrawEllipse(glm::vec3(dx, 0, 0), glm::vec3(00, e0.x, e0.y), glm::vec3(0, e1.x, e1.y), { 0, 128, 0 });
        }


        //static glm::vec3 P = { 0 ,0 ,0 };
        //float d = glm::dot( P, inv_cov * P );
        //ManipulatePosition(camera, &P, 0.2f);

        //char label[128];
        //sprintf(label, "%.2f", d);
        //DrawText(P, label);
        //

        for( float depth = -1.0f ; depth <= 1.0f ; depth += 0.03f )
        {
            float kPrime2 = 1.0f - sqr(glm::dot(W, W * depth) / glm::dot(W, W));
            float kPrime = std::sqrt(kPrime2);

            DrawEllipse(W * depth, U * kPrime, V * kPrime, { 128 , 128 , 128 }, 64);
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());

        ImGui::SliderFloat("fov", &camera.fovy, 0, 0.1);
        ImGui::Text("cam d %f", glm::length(camera.origin));

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
