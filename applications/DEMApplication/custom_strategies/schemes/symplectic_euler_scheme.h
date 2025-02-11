//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_SYMPLECTIC_EULER_SCHEME_H_INCLUDED )
#define  KRATOS_SYMPLECTIC_EULER_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes
#include "dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "custom_utilities/GeometryFunctions.h"
#include "utilities/quaternion.h"
#include "includes/ublas_interface.h"

namespace Kratos {

    namespace ublas = boost::numeric::ublas;

    class KRATOS_API(DEM_APPLICATION) SymplecticEulerScheme : public DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of SymplecticEulerScheme
        KRATOS_CLASS_POINTER_DEFINITION(SymplecticEulerScheme);

        /// Default constructor.
        SymplecticEulerScheme() {}

        /// Destructor.
        virtual ~SymplecticEulerScheme() {}

        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new SymplecticEulerScheme(*this));
            return cloned_scheme;
        }

        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new SymplecticEulerScheme(*this));
            return cloned_scheme;
        }

        // Function to calculate the inverse of a matrix using LU decomposition
        template<typename T>
        bool invert_matrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
            ublas::matrix<T> A(input);
            ublas::permutation_matrix<std::size_t> pm(A.size1());

            int res = ublas::lu_factorize(A, pm);
            if (res != 0) return false;

            inverse.assign(ublas::identity_matrix<T>(A.size1()));
            ublas::lu_substitute(A, pm, inverse);
            return true;
        }

        void SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override;
        void SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        void UpdateTranslationalVariables(
                int StepFlag,
                Node& i,
                array_1d<double, 3 >& coor,
                array_1d<double, 3 >& displ,
                array_1d<double, 3 >& delta_displ,
                array_1d<double, 3 >& vel,
                const array_1d<double, 3 >& initial_coor,
                const array_1d<double, 3 >& force,
                const double force_reduction_factor,
                const double mass,
                const double delta_t,
                const bool Fix_vel[3]) override;

        void CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        virtual void CalculateNewRotationalVariablesOfPolyhedrons(
                int StepFlag,
                Node& i,
                Matrix moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]);

        void CalculateNewRotationalVariablesOfRigidBodyElements(
                int StepFlag,
                Node& i,
                const array_1d<double, 3 > moments_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        void UpdateRotationalVariables(
                int StepFlag,
                Node& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) override;

        void CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) override;

        virtual void CalculateAngularAccelerationForPolyhedron(
                Matrix moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration);

        void CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration) override;

        /// Turn back information as a string.

        virtual std::string Info() const override {
            std::stringstream buffer;
            buffer << "SymplecticEulerScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override {
            rOStream << "SymplecticEulerScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override {
        }


    protected:


    private:

    /// Assignment operator.

        SymplecticEulerScheme& operator=(SymplecticEulerScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        SymplecticEulerScheme(SymplecticEulerScheme const& rOther) {
            *this = rOther;
        }

        ///@}

    }; // Class SymplecticEulerScheme

    inline std::istream& operator>>(std::istream& rIStream,
            SymplecticEulerScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const SymplecticEulerScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_SYMPLECTIC_EULER_SCHEME_H_INCLUDED  defined
