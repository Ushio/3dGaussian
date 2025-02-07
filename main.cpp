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

inline float ss_sqrt(float x) {
    float y; 
    _mm_store_ss(&y, _mm_sqrt_ss(_mm_load_ss(&x))); 
    return y; 
}

float min_of_abs(float x, float y, float z)
{
    return ss_min(ss_min(glm::abs(x), glm::abs(y)), glm::abs(z));
}
float max_of_abs(float x, float y)
{
    return ss_max(glm::abs(x), glm::abs(y));
}
void eigen_decomposition( glm::vec3 es[3], float lambdas[3], float A_00, float A_01, float A_02, float A_11, float A_12, float A_22 )
{
    glm::vec3 wbasisXY[2] = { {1, 0, 0}, {0, 1, 0} };

    auto sincostan = []( float* s, float* c, float* t, float cos_2theta, float sin_2theta )
    {                
        float X = cos_2theta;
        float Y = sin_2theta;
        float YY = Y * Y;
        float L = ss_sqrt(ss_max(X * X + YY, FLT_MIN) /* handle the singularity */);

        // The half vector
        float hx = X + sign_of(X) * L;
        float hy = Y;
        float hyhy = YY;
        float lh = ss_sqrt(hx * hx + hyhy);

        float tanTheta = hy / hx;
        float cosTheta = hx / lh;
        *t = tanTheta;
        *c = cosTheta;
        *s = cosTheta * tanTheta;
    };

    for (int i = 0; i < 32; i++)
    {
        {
            float b = A_12;
            float a = A_11;
            float d = A_22;

            // b is small enough not to affect the next rotation
            if( A_02 + b != A_02 )
            {
                float c;
                float s;
                float t;
                sincostan( &s, &c, &t, 0.5f * (d - a), b );

                //glm::mat3 P = glm::mat3(
                //    1, 0, 0,
                //    0, c, -s,
                //    0, s, c
                //);

                for (int j = 0; j < 2; j++)
                {
                    float Y = wbasisXY[j].y;
                    float Z = wbasisXY[j].z;
                    wbasisXY[j].y = c * Y - s * Z;
                    wbasisXY[j].z = s * Y + c * Z;
                }

                float nA_01 = c * A_01 /* - s * A_02*/;
                float nA_02 = /* c * A_02 */ + s * A_01;

                if( i == 0 )
                {
                    nA_01 += -s * A_02;
                    nA_02 += c * A_02;
                }

                // simplified via nA_12 == 0
                float nA_11 = A_11 - t * A_12;
                float nA_22 = A_22 + t * A_12;

                float nA_12 = 0.0f; // focus
                A_01 = nA_01;
                A_02 = nA_02;
                A_11 = nA_11;
                A_12 = nA_12;
                A_22 = nA_22;


                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = min_of_abs(A_00, A_11, A_22);
                if (minDiag + max_of_abs(A_02, A_01) == minDiag)
                {
                    break;
                }
            }
        }

        {
            float b = A_01;
            float a = A_00;
            float d = A_11;

            // b is small enough not to affect the next rotation
            if (A_12 + b != A_12)
            {
                float c;
                float s;
                float t;
                sincostan(&s, &c, &t, 0.5f * (d - a), b);

                //glm::mat3 P = glm::mat3(
                //    c, -s, 0,
                //    s, c, 0,
                //    0, 0, 1
                //);

                for (int j = 0; j < 2; j++)
                {
                    float X = wbasisXY[j].x;
                    float Y = wbasisXY[j].y;
                    wbasisXY[j].x = c * X - s * Y;
                    wbasisXY[j].y = s * X + c * Y;
                }

                // simplified via nA_12 == 0
                float nA_00 = A_00 - t * A_01;
                float nA_11 = A_11 + t * A_01;

                float nA_01 = 0.0f; // focus
                float nA_02 = c * A_02 /*- s * A_12*/;
                float nA_12 = /* c * A_12 + */ s * A_02;
                A_00 = nA_00;
                A_01 = nA_01;
                A_02 = nA_02;
                A_11 = nA_11;
                A_12 = nA_12;

                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = min_of_abs(A_00, A_11, A_22);
                if (minDiag + max_of_abs(A_02, A_12) == minDiag)
                {
                    break;
                }
            }
        }

        {
            float b = A_02;
            float a = A_00;
            float d = A_22;
            // b is small enough not to affect the next rotation
            if( A_01 + b != A_01 )
            {
                float c;
                float s;
                float t;
                sincostan(&s, &c, &t, 0.5f * (d - a), b);

                //glm::mat3 P = glm::mat3(
                //    c, 0, -s,
                //    0, 1, 0,
                //    s, 0, c
                //);

                for (int j = 0; j < 2; j++)
                {
                    float X = wbasisXY[j].x;
                    float Z = wbasisXY[j].z;
                    wbasisXY[j].x = c * X - s * Z;
                    wbasisXY[j].z = s * X + c * Z;
                }

                // simplified via nA_12 == 0
                float nA_00 = A_00 - t * A_02;
                float nA_22 = A_22 + t * A_02;

                float nA_01 = /* c * A_01 */ - s * A_12;
                float nA_02 = 0.0f; // focus
                float nA_12 = c * A_12 /* + s * A_01*/;
                A_00 = nA_00;
                A_01 = nA_01;
                A_02 = nA_02;
                A_12 = nA_12;
                A_22 = nA_22;


                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = min_of_abs(A_00, A_11, A_22);
                if (minDiag + max_of_abs(A_12, A_01) == minDiag)
                {
                    break;
                }
            }
        }
    }

    glm::vec3 wbasisZ = glm::cross(wbasisXY[0], wbasisXY[1]);
    lambdas[0] = A_00;
    lambdas[1] = A_11;
    lambdas[2] = A_22;
    es[0] = { wbasisXY[0].x, wbasisXY[1].x, wbasisZ.x };
    es[1] = { wbasisXY[0].y, wbasisXY[1].y, wbasisZ.y };
    es[2] = { wbasisXY[0].z, wbasisXY[1].z, wbasisZ.z };
}

// single side jacobi compact implementation Assume B is symmetric
void eigen_decomposition_jacobi(glm::vec3 es[3], float lambdas[3], glm::mat3 B )
{
    float convergence_previous = FLT_MAX;

    for(;;)
    {
        float convergence = 0.0f;

        for (int i = 0; i < 3; i++)
        {
            int index_b1 = i;
            int index_b2 = (i + 1) % 3;
            glm::vec3 b1 = B[index_b1];
            glm::vec3 b2 = B[index_b2];

            float Py = 2.0f * glm::dot(b1, b2);
            float Px = glm::dot(b1, b1) - glm::dot(b2, b2);
            float PL = sqrtf(Px * Px + Py * Py);

            convergence = glm::max(convergence, glm::abs(Py));

            if (PL == 0.0f)
            {
                continue; // no rotation
            }

            float sgn = 0.0f < Px ? 1.0f : -1.0f;
            glm::vec2 H = { PL + Px * sgn, Py * sgn };
            float L = glm::length(H);
            float c = H.x / L;
            float s = H.y / L;

            // Equivalent to:
            //float two_theta = atan(Py / Px);
            //float s = sin(two_theta * 0.5f);
            //float c = cos(two_theta * 0.5f);

            for (int j = 0; j < 3; j++)
            {
                B[index_b1][j] = +c * b1[j] + s * b2[j];
                B[index_b2][j] = -s * b1[j] + c * b2[j];
            }
        }

        if (convergence < convergence_previous && convergence != 0.0f )
        {
            convergence_previous = convergence;
            continue;
        }
        break;
    }

    for (int i = 0; i < 3; i++)
    {
        float l = glm::length(B[i]);
        lambdas[i] = l;
        es[i] = B[i] / l;
    }
}

glm::vec3 plasma_quintic(float x)
{
    x = glm::clamp(x, 0.0f, 1.0f);
    glm::vec4 x1 = glm::vec4(1.0f, x, x * x, x * x * x); // 1 x x2 x3
    glm::vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
    return glm::vec3(
        glm::dot(x1, glm::vec4(+0.063861086f, +1.992659096f, -1.023901152f, -0.490832805f)) + glm::dot(glm::vec2(x2), glm::vec2(+1.308442123f, -0.914547012f)),
        glm::dot(x1, glm::vec4(+0.049718590f, -0.791144343f, +2.892305078f, +0.811726816f)) + glm::dot(glm::vec2(x2), glm::vec2(-4.686502417f, +2.717794514f)),
        glm::dot(x1, glm::vec4(+0.513275779f, +1.580255060f, -5.164414457f, +4.559573646f)) + glm::dot(glm::vec2(x2), glm::vec2(-1.916810682f, +0.570638854f)));
}

glm::uvec4 pcg4d(glm::uvec4 v) {
    v = v * 1664525u + 1013904223u; 
    v.x += v.y * v.w;
    v.y += v.z * v.x;
    v.z += v.x * v.y;
    v.w += v.y * v.z;
    v ^= v >> 16u;
    v.x += v.y * v.w;
    v.y += v.z * v.x;
    v.z += v.x * v.y;
    v.w += v.y * v.z;
    return v; 
}
float uniformf(uint32_t x)
{
    uint32_t bits = (x >> 9) | 0x3f800000;
    float value;
    memcpy(&value, &bits, sizeof(float));
    return value - 1.0f;
}

glm::vec4 boxmular4( float X0, float Y0, float X1, float Y1 )
{
    float k0 = ss_min( sqrt(-2.0f * log(X0)), FLT_MAX );
    float k1 = ss_min( sqrt(-2.0f * log(X1)), FLT_MAX );
    float Z0 = k0 * cos(2.0f * glm::pi<float>() * Y0);
    float Z1 = k0 * sin(2.0f * glm::pi<float>() * Y0);
    float Z2 = k1 * cos(2.0f * glm::pi<float>() * Y1);
    float Z3 = k1 * sin(2.0f * glm::pi<float>() * Y1);
    return { Z0, Z1, Z2, Z3 };
}

int main() {
    using namespace pr;

    // simple tests
    Xoshiro128StarStar random;
    Stopwatch sw;
    for (int i = 0; i < 1000000; i++)
    {
        glm::vec3 U = {
            glm::mix( -100.0f, 100.0f, random.uniformf() ),
            glm::mix( -100.0f, 100.0f, random.uniformf() ),
            glm::mix( -100.0f, 100.0f, random.uniformf() ),
        };
        glm::vec3 V = {
            glm::mix(-100.0f, 100.0f, random.uniformf()),
            glm::mix(-100.0f, 100.0f, random.uniformf()),
            glm::mix(-100.0f, 100.0f, random.uniformf()),
        };
        glm::vec3 W = {
            glm::mix(-100.0f, 100.0f, random.uniformf()),
            glm::mix(-100.0f, 100.0f, random.uniformf()),
            glm::mix(-100.0f, 100.0f, random.uniformf()),
        };

        {
            glm::vec3 N = glm::cross(U, V);
            V = glm::normalize(glm::cross(N, U)) * glm::length(V);
            W = glm::normalize(N) * glm::length(W);
        }

        glm::mat3 R = glm::mat3(glm::normalize(U), glm::normalize(V), glm::normalize(W));
        glm::mat3 L = glm::mat3(
            lengthSquared(U), 0, 0,
            0, lengthSquared(V), 0,
            0, 0, lengthSquared(W)
        );
        glm::mat3 cov = R * L * glm::transpose(R);

        float lambdas[3];
        glm::vec3 es[3];
        //eigen_decomposition(es, lambdas, cov[0][0], cov[1][0], cov[2][0], cov[1][1], cov[2][1], cov[2][2]);

        eigen_decomposition_jacobi(es, lambdas, cov);

        //printf("%f\n", glm::abs(glm::dot(es[1], es[2])));
        PR_ASSERT(glm::abs(glm::dot(es[0], es[1])) < 1.0e-6f);
        PR_ASSERT(glm::abs(glm::dot(es[1], es[2])) < 1.0e-6f);
        PR_ASSERT(glm::abs(glm::dot(es[2], es[0])) < 1.0e-6f);

        float S_eigen_ref = ss_min(ss_min(L[0][0], L[1][1]), L[2][2]);
        float L_eigen_ref = ss_max(ss_max(L[0][0], L[1][1]), L[2][2]);
        float S_eigen = ss_min(ss_min(lambdas[0], lambdas[1]), lambdas[2]);
        float L_eigen = ss_max(ss_max(lambdas[0], lambdas[1]), lambdas[2]);
        float epsL = FLT_EPSILON * L_eigen_ref * 6.0f;
        PR_ASSERT( glm::abs( S_eigen - S_eigen_ref ) < epsL);
        PR_ASSERT( glm::abs( L_eigen - L_eigen_ref ) < epsL);

        float meanTr = (cov[0][0] + cov[1][1] + cov[2][2]) / 3.0f;
        float meanLambda = (lambdas[0] + lambdas[1] + lambdas[2]) / 3.0f;
        PR_ASSERT(glm::abs(meanTr - meanLambda) < FLT_EPSILON * meanLambda * 8.0f);
    }
    printf("%f\n", sw.elapsed() * 1000.0f);
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
    // camera.fovy = glm::radians( 90.0f );

    //camera.origin *= 100.0f;
    //camera.fovy = 0.005f;

    double e = GetElapsedTime();

    ITexture* tex = CreateTexture();
    Image2DRGBA8 image;
    int stride = 1;

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        // ClearBackground(0.1f, 0.1f, 0.1f, 1);
        ClearBackground(tex);

        BeginCamera(camera);

        PushGraphicState();

        DrawGrid(GridAxis::XZ, 1.0f, 10, { 128, 128, 128 });

        static glm::vec3 U = { 1, 0, 0 };
        static glm::vec3 V = { 0, 1, 0 };
        static glm::vec3 W = { 0, 0, 1 };
        //static glm::vec3 U = { 0.363252401, 0.425909996, -0.522605419 };
        //static glm::vec3 V = { -0.283039778, 0.830360413, 0.830360413 };
        //static glm::vec3 W = { 0.833595634, -0.0345231034,0.551279962 };
        ManipulatePosition(camera, &U, 0.2f);
        ManipulatePosition(camera, &V, 0.2f);
        ManipulatePosition(camera, &W, 0.2f);

        static bool RT_mode = false;
        static bool showWire = true;

        //{
        //    glm::vec3 N = glm::cross(U, V);
        //    V = glm::normalize(glm::cross(N, U)) * glm::length(V);

        //    DrawArrow({ 0,0,0 }, U, 0.01f, { 255,0,0 });
        //    DrawArrow({ 0,0,0 }, V, 0.01f, { 0,255,0 });
        //    DrawArrow({ 0,0,0 }, W, 0.01f, { 0,0,255 });

        //    W = glm::normalize(N) * glm::length(W);
        //}

        {
            float cov_00 = U.x * U.x + V.x * V.x + W.x * W.x;
            float cov_01 = U.x * U.y + V.x * V.y + W.x * W.y;
            float cov_02 = U.x * U.z + V.x * V.z + W.x * W.z;
            float cov_11 = U.y * U.y + V.y * V.y + W.y * W.y;
            float cov_12 = U.y * U.z + V.y * V.z + W.y * W.z;
            float cov_22 = U.z * U.z + V.z * V.z + W.z * W.z;
            glm::mat3 cov = {
                cov_00, cov_01, cov_02,
                cov_01, cov_11, cov_12,
                cov_02, cov_12, cov_22
            };

            // https://tavianator.com/2014/ellipsoid_bounding_boxes.html
            float dx = std::sqrt(cov_00);// std::sqrt(U.x * U.x + V.x * V.x + W.x * W.x);
            float dy = std::sqrt(cov_11);// std::sqrt(U.y * U.y + V.y * V.y + W.y * W.y);
            float dz = std::sqrt(cov_22);// std::sqrt(U.z * U.z + V.z * V.z + W.z * W.z);
            //DrawCube({ 0, 0, 0 }, { dx * 2, dy * 2, dz * 2 }, { 255,255,255 });

            float lambdas[3];
            glm::vec3 es[3];
            // eigen_decomposition(es, lambdas, cov[0][0], cov[1][0], cov[2][0], cov[1][1], cov[2][1], cov[2][2]);
            eigen_decomposition_jacobi(es, lambdas, cov);


            DrawArrow({ 0,0,0 }, es[0], 0.01f, { 255,0,255 });
            DrawArrow({ 0,0,0 }, es[1], 0.01f, { 255,255,0 });
            DrawArrow({ 0,0,0 }, es[2], 0.01f, { 0,255,255 });
            DrawLine(-es[0] * std::sqrt(lambdas[0]), es[0] * std::sqrt(lambdas[0]), { 255,0,255 });
            DrawLine(-es[1] * std::sqrt(lambdas[1]), es[1] * std::sqrt(lambdas[1]), { 255,255,0 });
            DrawLine(-es[2] * std::sqrt(lambdas[2]), es[2] * std::sqrt(lambdas[2]), { 0,255,255 });

            if( showWire )
            {
                // glm::mat3 ellipsoidAffine = { U, V, W };
                //SetObjectTransform(ellipsoidAffine);
                //DrawSphere({ 0,0,0 }, 1.0f, { 128 ,128 ,128 }, 32, 32);
                //SetObjectIdentify();
                glm::vec3 e0 = es[0] * std::sqrt(lambdas[0]);
                glm::vec3 e1 = es[1] * std::sqrt(lambdas[1]);
                glm::vec3 e2 = es[2] * std::sqrt(lambdas[2]);
                DrawEllipse({ 0,0,0 }, e0, e1, { 255,255,255 }, 64);
                DrawEllipse({ 0,0,0 }, e1, e2, { 255,255,255 }, 64);
                DrawEllipse({ 0,0,0 }, e2, e0, { 255,255,255 }, 64);
            }
        }

        //glm::mat3 ellipsoidAffine = { U, V, W };

        //for (int i = 0; i < 10000; i++)
        //{
        //    glm::uvec4 r = pcg4d({ i, 0, 0, 0 });
        //    glm::vec3 p = boxmular4(uniformf(r.x), uniformf(r.y), uniformf(r.z), uniformf(r.w));

        //    p = ellipsoidAffine * p;

        //    DrawSphere(p, 0.01f, { 255, 255, 255 }, 4, 4);
        //}


        image.allocate(GetScreenWidth() / stride, GetScreenHeight() / stride);
        CameraRayGenerator rayGenerator(GetCurrentViewMatrix(), GetCurrentProjMatrix(), image.width(), image.height());

        glm::mat4 viewMat = GetCurrentViewMatrix();
        glm::mat3 viewRot = glm::mat3(viewMat);

        {
            glm::vec3 u_ray = viewMat * glm::vec4(0, 0, 0, 1);
            float zcomp0 = -u_ray.x / (u_ray.z * u_ray.z);
            float zcomp1 = -u_ray.y / (u_ray.z * u_ray.z);
            printf("%f %f\n", zcomp0, zcomp1);
        }

        if(0)
        for (int j = 0; j < image.height(); ++j)
        {
            for (int i = 0; i < image.width(); ++i)
            {
                if (RT_mode)
                {
                    glm::vec3 ro, rd;
                    rayGenerator.shoot(&ro, &rd, i, j, 0.5f, 0.5f);

                    glm::vec3 Z = glm::normalize(rd);
                    glm::vec3 X;
                    glm::vec3 Y;
                    GetOrthonormalBasis(Z, &X, &Y);

                    glm::mat3 toLocal = {
                        glm::vec3(X.x, Y.x, Z.x),
                        glm::vec3(X.y, Y.y, Z.y),
                        glm::vec3(X.z, Y.z, Z.z),
                    };
                    glm::vec3 u = toLocal * U;
                    glm::vec3 v = toLocal * V;
                    glm::vec3 w = toLocal * W;

                    glm::vec2 v_rel = toLocal * ro;

                    float cov_00 = u.x * u.x + v.x * v.x + w.x * w.x;
                    float cov_01 = u.x * u.y + v.x * v.y + w.x * w.y;
                    float cov_11 = u.y * u.y + v.y * v.y + w.y * w.y;
                    glm::mat2 cov = { cov_00, cov_01, cov_01, cov_11 };
                    glm::mat2 inv_cov = glm::inverse(cov);

                    float d2 = glm::dot(v_rel, inv_cov * v_rel);
                    float alpha = glm::exp(-0.5f * d2);
                    glm::u8vec3 color = glm::u8vec3(glm::clamp(plasma_quintic(alpha) * 255.0f, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(255.0f, 255.0f, 255.0f)));
                
                    if (glm::abs(d2 - 1.0f) < 0.005f)
                    {
                        color = { 255,255,255 };
                    }
                    if (glm::abs(d2 - 4.0f) < 0.005f)
                    {
                        color = { 128,128,128 };
                    }
                    image(i, j) = glm::u8vec4(color, 255);
                }
                else
                {
                    // viewRot
                    glm::vec3 u = viewRot * U;
                    glm::vec3 v = viewRot * V;
                    glm::vec3 w = viewRot * W;

                    float cov_00 = u.x * u.x + v.x * v.x + w.x * w.x;
                    float cov_01 = u.x * u.y + v.x * v.y + w.x * w.y;
                    float cov_02 = u.x * u.z + v.x * v.z + w.x * w.z;
                    float cov_11 = u.y * u.y + v.y * v.y + w.y * w.y;
                    float cov_12 = u.y * u.z + v.y * v.z + w.y * w.z;
                    float cov_22 = u.z * u.z + v.z * v.z + w.z * w.z;
                    glm::mat3 cov = {
                        cov_00, cov_01, cov_02,
                        cov_01, cov_11, cov_12,
                        cov_02, cov_12, cov_22
                    };
                    // glm::vec3 u_ray = viewMat * glm::vec4( ( /* gaussian center */ - camera.origin ), 1.0f );
                    glm::vec3 u_ray = viewMat * glm::vec4( 0, 0, 0, 1 );

                    glm::mat3 J;
                    J = {
                        1.0f / u_ray.z, 0, 0,
                        0, 1.0f / u_ray.z, 0,
                        -u_ray.x / (u_ray.z * u_ray.z), -u_ray.y / (u_ray.z * u_ray.z), 1.0f
                    };

                    glm::mat3 covPrime = J * cov * glm::transpose( J );
                    glm::mat2 covPrime2d = glm::mat2(
                        covPrime[0][0], covPrime[0][1],
                        covPrime[1][0], covPrime[1][1]
                    );

                    glm::mat2 inv_cov = glm::inverse(covPrime2d);

                    float tanThetaY = std::tan( camera.fovy * 0.5f );
                    float tanThetaX = tanThetaY / image.height() * image.width();

                    glm::vec2 v_rel = glm::vec2(
                        glm::mix(-tanThetaX, tanThetaX, (float)i / image.width()),
                        glm::mix(tanThetaY, -tanThetaY, (float)j / image.height())
                    ) - glm::vec2(u_ray.x / -u_ray.z, u_ray.y / -u_ray.z); // z is negative

                    float d2 = glm::dot(v_rel, inv_cov * v_rel);
                    float alpha = std::expf(-0.5f * d2);
                    glm::u8vec3 color = glm::u8vec3(glm::clamp(plasma_quintic(alpha) * 255.0f, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(255.0f, 255.0f, 255.0f)));

                    if (glm::abs(d2 - 1.0f) < 0.005f)
                    {
                        color = { 255,255,255 };
                    }
                    if (glm::abs(d2 - 4.0f) < 0.005f)
                    {
                        color = { 128,128,128 };
                    }
                    image(i, j) = glm::u8vec4(color, 255);

                    // image(i, j) = glm::u8vec4(0, 0, 0, 255);
                }
            }
        }
        tex->upload(image);

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());

        ImGui::SliderFloat("fov", &camera.fovy, 0, glm::radians( 170.0f));
        ImGui::Text("cam d %f", glm::length(camera.origin));
        ImGui::Checkbox("RT_mode", &RT_mode);
        ImGui::Checkbox("showWire", &showWire);
        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
