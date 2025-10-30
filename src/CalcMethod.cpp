#include "CalcMethod.h"
#include "utils/ArrayUtils.h"

void StormerVerletMethod::calculateX(const double dt) {
    for (auto& p : particles) {
        const auto x_curr = p.getX();
        const double m = p.getM();
        const auto F = p.getF();
        const auto v = p.getV();
        const auto a = (1.0 / m) * F;
        const auto x_new = x_curr + dt * v + 0.5 * (dt * dt) * a;
        p.setX(x_new);
    }
}

void StormerVerletMethod::calculateV(const double dt) {
    for (auto& p : particles) {
        const auto v_curr = p.getV();
        const double m_i = p.getM();
        const auto F = p.getF();
        const auto F_old = p.getOldF();
        const auto v_new = v_curr + (dt / (2.0 * m_i)) * (F_old + F);
        p.setV(v_new);
    }
}

void StormerVerletMethod::calculateF(const double dt) {
    for (auto& p : particles) {
        p.setF({});
    }

    const size_t n_particles = particles.size();
    for (size_t i = 0; i < n_particles; ++i) {
        // index offset for Newton's third law
        for (size_t j = i + 1; j < n_particles; ++j) {
            constexpr double norm_squared_soft_const = 1e-9;
            auto& p_i = particles[i];
            auto& p_j = particles[j];

            const auto dist = p_j.getX() - p_i.getX();
            const double norm = ArrayUtils::L2Norm(dist);
            if (norm < EPSILON) {
                // avoid division by zero and numerical explosion when particles are at the same position or very close
                continue;
            }
            const double norm_squared_softened =
                norm * norm + norm_squared_soft_const;
            const double norm_cubed_softened = norm_squared_softened * norm;

            const auto F_vector =
                ((p_i.getM() * p_j.getM()) / norm_cubed_softened) * dist;

            // apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
            auto F_i = p_i.getF();
            auto F_j = p_j.getF();
            // actio est reactio
            F_i = F_i + F_vector;
            F_j = F_j - F_vector;
            p_i.setF(F_i);
            p_j.setF(F_j);
        }
    }
}
