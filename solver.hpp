/*
 * solver.hpp
 *
 *  Created on: Jan 3, 2021
 *      Author: d-w-h
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include "user_types.hpp"

void solver(d_data domain_data,
            t_data time_data,
            p_params physical_params,
            b_data boundaries_and_psi_init,
            s_data* solver_data);

#endif /* SOLVER_HPP_ */
