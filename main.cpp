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

void eigen_decompress( glm::vec3 es[3], float lambdas[3], glm::mat3 M )
{
    float A_00 = M[0][0];
    float A_01 = M[1][0];
    float A_02 = M[2][0];
    float A_11 = M[1][1];
    float A_12 = M[2][1];
    float A_22 = M[2][2];

    glm::vec3 wbasisXY[2] = { {1, 0, 0}, {0, 1, 0} };

    auto sincos_of = [](float* s, float* c, float invTan2Theta)
    {
        float tanTheta = 1.0f / (sign_of(invTan2Theta) * std::sqrtf(1.0f + invTan2Theta * invTan2Theta) + invTan2Theta);
        float cosTheta = 1.0f / std::sqrtf(1.0f + tanTheta * tanTheta);
        float sinTheta = tanTheta * cosTheta;
        *s = sinTheta;
        *c = cosTheta;
    };
    const float eps_nondiag = 1.0e-15f;
    for (int i = 0; i < 32; i++)
    {
        {
            float b = A_12;
            float a = A_11;
            float d = A_22;

            if( eps_nondiag < glm::abs(b) )
            {
                float c;
                float s;
                float invTan2Theta = 0.5f * (d - a) / b;
                sincos_of(&s, &c, invTan2Theta);

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

                float nA_11 = c * (c * A_11 - s * A_12) - s * (c * A_12 - s * A_22);
                float nA_12 = 0.0f; // focus
                float nA_22 = s * (c * A_12 + s * A_11) + c * (c * A_22 + s * A_12);
                A_01 = nA_01;
                A_02 = nA_02;
                A_11 = nA_11;
                A_12 = nA_12;
                A_22 = nA_22;

                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = ss_min(ss_min(glm::abs(A_00), glm::abs(A_11)), glm::abs(A_22));
                if (minDiag + ss_max(glm::abs(A_02), glm::abs(A_01)) == minDiag)
                {
                    break;
                }
            }
        }

        {
            float b = A_01;
            float a = A_00;
            float d = A_11;
            if (eps_nondiag < glm::abs(b))
            {
                float c;
                float s;
                float invTan2Theta = 0.5f * (d - a) / b;
                sincos_of(&s, &c, invTan2Theta);

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

                float nA_00 = c * (c * A_00 - s * A_01) - s * (c * A_01 - s * A_11);
                float nA_01 = 0.0f; // focus
                float nA_02 = c * A_02 /*- s * A_12*/;
                float nA_11 = s * (c * A_01 + s * A_00) + c * (c * A_11 + s * A_01);
                float nA_12 = /* c * A_12 + */ s * A_02;
                A_00 = nA_00;
                A_01 = nA_01;
                A_02 = nA_02;
                A_11 = nA_11;
                A_12 = nA_12;

                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = ss_min(ss_min(glm::abs(A_00), glm::abs(A_11)), glm::abs(A_22));
                if (minDiag + ss_max(glm::abs(A_02), glm::abs(A_12)) == minDiag)
                {
                    break;
                }
            }
        }

        {
            float b = A_02;
            float a = A_00;
            float d = A_22;
            if (eps_nondiag < glm::abs(b))
            {
                float c;
                float s;
                float invTan2Theta = 0.5f * (d - a) / b;
                sincos_of(&s, &c, invTan2Theta);

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

                float nA_00 = c * (c * A_00 - s * A_02) - s * (c * A_02 - s * A_22);
                float nA_01 = /* c * A_01 */ - s * A_12;
                float nA_02 = 0.0f; // focus
                float nA_12 = c * A_12 /* + s * A_01*/;
                float nA_22 = s * (c * A_02 + s * A_00) + c * (c * A_22 + s * A_02);
                A_00 = nA_00;
                A_01 = nA_01;
                A_02 = nA_02;
                A_12 = nA_12;
                A_22 = nA_22;

                // Converged when the smallest diag can not get affected from the largest non-diag
                float minDiag = ss_min(ss_min(glm::abs(A_00), glm::abs(A_11)), glm::abs(A_22));
                if (minDiag + ss_max(glm::abs(A_12), glm::abs(A_01)) == minDiag)
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

int main() {
    using namespace pr;

    // simple tests
    Xoshiro128StarStar random;
    for (int i = 0; i < 10000; i++)
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
        eigen_decompress(es, lambdas, cov);

        PR_ASSERT(glm::abs(glm::dot(es[0], es[1])) < 1.0e-6f);
        PR_ASSERT(glm::abs(glm::dot(es[1], es[2])) < 1.0e-6f);
        PR_ASSERT(glm::abs(glm::dot(es[2], es[0])) < 1.0e-6f);

        float S_eigen_ref = ss_min(ss_min(L[0][0], L[1][1]), L[2][2]);
        float L_eigen_ref = ss_max(ss_max(L[0][0], L[1][1]), L[2][2]);
        float S_eigen = ss_min(ss_min(lambdas[0], lambdas[1]), lambdas[2]);
        float L_eigen = ss_max(ss_max(lambdas[0], lambdas[1]), lambdas[2]);
        float eps0 = FLT_EPSILON * S_eigen_ref * 256.0f;
        float eps1 = FLT_EPSILON * L_eigen_ref * 256.0f;
        PR_ASSERT( glm::abs( S_eigen - S_eigen_ref ) < eps0);
        PR_ASSERT( glm::abs( L_eigen - L_eigen_ref ) < eps1);

        float det = glm::determinant(cov);
        PR_ASSERT(glm::abs(lambdas[0] * lambdas[1] * lambdas[2] - det ) < FLT_EPSILON * det * 256.0f);

        float meanTr = (cov[0][0] + cov[1][1] + cov[2][2]) / 3.0f;
        float meanLambda = (lambdas[0] + lambdas[1] + lambdas[2]) / 3.0f;
        PR_ASSERT(glm::abs(meanTr - meanLambda) < FLT_EPSILON * meanLambda * 8.0f);
    }

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

        //static glm::vec3 U = { 1, 0, 0 };
        //static glm::vec3 V = { 0, 1, 0 };
        //static glm::vec3 W = { 0, 0, 1 };
        static glm::vec3 U = { 0.363252401, 0.425909996, -0.522605419 };
        static glm::vec3 V = { -0.283039778, 0.830360413, 0.830360413 };
        static glm::vec3 W = { 0.833595634, -0.0345231034,0.551279962 };
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

        float lambdas[3];
        glm::vec3 es[3];
        eigen_decompress(es, lambdas, cov);
        DrawArrow({ 0,0,0 }, es[0], 0.01f, { 255,0,255 });
        DrawArrow({ 0,0,0 }, es[1], 0.01f, { 255,255,0 });
        DrawArrow({ 0,0,0 }, es[2], 0.01f, { 0,255,255 });
        DrawLine(-es[0] * std::sqrt(lambdas[0]), es[0] * std::sqrt(lambdas[0]), { 255,0,255 });
        DrawLine(-es[1] * std::sqrt(lambdas[1]), es[1] * std::sqrt(lambdas[1]), { 255,255,0 });
        DrawLine(-es[2] * std::sqrt(lambdas[2]), es[2] * std::sqrt(lambdas[2]), { 0,255,255 });

        // naiive 
        //{
        //    //glm::mat3 V_all = glm::identity<glm::mat3>();
        //    float A_00 = cov[0][0];
        //    float A_01 = cov[1][0];
        //    float A_02 = cov[2][0];
        //    float A_11 = cov[1][1];
        //    float A_12 = cov[2][1];
        //    float A_22 = cov[2][2];

        //    glm::vec3 wbasisXY[2] = { {1, 0, 0}, {0, 1, 0} };

        //    auto sincos_of = []( float *s, float *c, float invTan2Theta )
        //    {
        //        float tanTheta = 1.0f / (sign_of(invTan2Theta) * std::sqrtf(1.0f + invTan2Theta * invTan2Theta) + invTan2Theta);
        //        float cosTheta = 1.0f / std::sqrtf(1.0f + tanTheta * tanTheta);
        //        float sinTheta = tanTheta * cosTheta;
        //        *s = sinTheta;
        //        *c = cosTheta;
        //    };

        //    const float eps = 1.0e-15f;
        //    for (int i = 0; i < nJacobiItr; i++)
        //    {
        //        {
        //            float b = A_12;
        //            float a = A_11;
        //            float d = A_22;

        //            if( eps < glm::abs( b ) )
        //            {
        //                float c;
        //                float s;
        //                float invTan2Theta = 0.5f * ( d - a ) / b;
        //                sincos_of( &s, &c, invTan2Theta );

        //                //glm::mat3 P = glm::mat3(
        //                //    1, 0, 0,
        //                //    0, c, -s,
        //                //    0, s, c
        //                //);

        //                for (int j = 0; j < 2; j++)
        //                {
        //                    float y = c * wbasisXY[j].y - s * wbasisXY[j].z;
        //                    float z = s * wbasisXY[j].y + c * wbasisXY[j].z;
        //                    wbasisXY[j].y = y;
        //                    wbasisXY[j].z = z;
        //                }

        //                float nA_01 = c * A_01 - s * A_02;
        //                float nA_02 = c * A_02 + s * A_01;
        //                float nA_11 = c * ( c * A_11 - s * A_12 ) - s * ( c * A_12 - s * A_22 );
        //                float nA_12 = 0.0f; // focus
        //                float nA_22 = s * ( c * A_12 + s * A_11 ) + c * ( c * A_22 + s * A_12 );
        //                A_01 = nA_01;
        //                A_02 = nA_02;
        //                A_11 = nA_11;
        //                A_12 = nA_12;
        //                A_22 = nA_22;
        //                //printf("%f\n", L[1][0] - nA_01);
        //                //printf("%f\n", L[2][0] - nA_02);
        //                //printf("%f\n", L[1][1] - nA_11);
        //                //printf("%f\n", L[2][2] - nA_22);
        //            }
        //        }

        //        {
        //            float b = A_01;
        //            float a = A_00;
        //            float d = A_11;
        //            if (eps < glm::abs(b))
        //            {
        //                float c;
        //                float s;
        //                float invTan2Theta = 0.5f * (d - a) / b;
        //                sincos_of(&s, &c, invTan2Theta);

        //                //glm::mat3 P = glm::mat3(
        //                //    c, -s, 0,
        //                //    s, c, 0,
        //                //    0, 0, 1
        //                //);

        //                for (int j = 0; j < 2; j++)
        //                {
        //                    float x = c * wbasisXY[j].x - s * wbasisXY[j].y;
        //                    float y = s * wbasisXY[j].x + c * wbasisXY[j].y;
        //                    wbasisXY[j].x = x;
        //                    wbasisXY[j].y = y;
        //                }

        //                float nA_00 = c * (c * A_00 - s * A_01) - s * (c * A_01 - s * A_11);
        //                float nA_01 = 0.0f; // focus
        //                float nA_02 = c * A_02 - s * A_12;
        //                float nA_11 = s * (c * A_01 + s * A_00) + c * (c * A_11 + s * A_01);
        //                float nA_12 = c * A_12 + s * A_02;
        //                A_00 = nA_00;
        //                A_01 = nA_01;
        //                A_02 = nA_02;
        //                A_11 = nA_11;
        //                A_12 = nA_12;
        //                //printf("%f\n", L[0][0] - nA_00);
        //                //printf("%f\n", L[1][0] - nA_01);
        //                //printf("%f\n", L[2][0] - nA_02);
        //                //printf("%f\n", L[1][1] - nA_11);
        //                //printf("%f\n", L[2][1] - nA_12);
        //            }
        //        }

        //        {
        //            float b = A_02;
        //            float a = A_00;
        //            float d = A_22;
        //            if (eps < glm::abs(b))
        //            {
        //                float c;
        //                float s;
        //                float invTan2Theta = 0.5f * (d - a) / b;
        //                sincos_of(&s, &c, invTan2Theta);

        //                //glm::mat3 P = glm::mat3(
        //                //    c, 0, -s,
        //                //    0, 1, 0,
        //                //    s, 0, c
        //                //);

        //                for (int j = 0; j < 2; j++)
        //                {
        //                    float x = c * wbasisXY[j].x - s * wbasisXY[j].z;
        //                    float z = s * wbasisXY[j].x + c * wbasisXY[j].z;
        //                    wbasisXY[j].x = x;
        //                    wbasisXY[j].z = z;
        //                }

        //                float nA_00 = c * ( c * A_00 - s * A_02 ) - s * ( c * A_02 - s * A_22 );
        //                float nA_01 = c * A_01 - s * A_12;
        //                float nA_02 = 0.0f; // focus
        //                float nA_12 = c * A_12 + s * A_01;
        //                float nA_22 = s * ( c * A_02 + s * A_00 ) + c * ( c * A_22 + s * A_02 );
        //                A_00 = nA_00;
        //                A_01 = nA_01;
        //                A_02 = nA_02;
        //                A_12 = nA_12;
        //                A_22 = nA_22;
        //                //printf("%f\n", L[0][0] - nA_00);
        //                //printf("%f\n", L[1][0] - nA_01);
        //                //printf("%f\n", L[2][0] - nA_02);
        //                //printf("%f\n", L[2][1] - nA_12);
        //                //printf("%f\n", L[2][2] - nA_22);
        //            }
        //        }
        //    }

        //    //DrawArrow({ 0,0,0 }, eigen0, 0.01f, { 255,0,255 });
        //    //DrawArrow({ 0,0,0 }, eigen1, 0.01f, { 255,255,0 });
        //    //DrawArrow({ 0,0,0 }, eigen2, 0.01f, { 0,255,255 });
        //    //DrawArrow({ 0,0,0 }, eigen0 * std::sqrt(A_00), 0.01f, { 255,0,255 });
        //    //DrawArrow({ 0,0,0 }, eigen1 * std::sqrt(A_11), 0.01f, { 255,255,0 });
        //    //DrawArrow({ 0,0,0 }, eigen2 * std::sqrt(A_22), 0.01f, { 0,255,255 });

        //    glm::vec3 wbasisZ = glm::cross(wbasisXY[0], wbasisXY[1]);
        //    glm::vec3 e0 = { wbasisXY[0].x, wbasisXY[1].x, wbasisZ.x };
        //    glm::vec3 e1 = { wbasisXY[0].y, wbasisXY[1].y, wbasisZ.y };
        //    glm::vec3 e2 = { wbasisXY[0].z, wbasisXY[1].z, wbasisZ.z };

        //    DrawArrow({ 0,0,0 }, e0, 0.01f, { 255,0,255 });
        //    DrawArrow({ 0,0,0 }, e1, 0.01f, { 255,255,0 });
        //    DrawArrow({ 0,0,0 }, e2, 0.01f, { 0,255,255 });
        //    DrawLine(-e0 * std::sqrt(A_00), e0 * std::sqrt(A_00), { 255,0,255 });
        //    DrawLine(-e1 * std::sqrt(A_11), e1 * std::sqrt(A_11), { 255,255,0 });
        //    DrawLine(-e2 * std::sqrt(A_22), e2 * std::sqrt(A_22), { 0,255,255 });
        //}

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
