// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_OPENCL_SAXS_DEBYE_H_
#define BCL_OPENCL_SAXS_DEBYE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_euclidean_distance.h"
#include "bcl_opencl_tools.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "restraint/bcl_restraint_sas_scattering_point.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SaxsDebye
    //! @brief Uses the Debye formula to calculate I from the form factors for a protein model
    //! @details I(q) = (sum i=1 to M)*(sum j=1 to M) Fi(q)*Fj(q) * sin(q*rij)/(q*rij)
    //!
    //! @see @link example_opencl_saxs_debye.cpp @endlink
    //! @author loweew
    //! @date May 5, 2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SaxsDebye :
      public restraint::SasDebyeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! bool value to represent loops that are not present in the protein model
      bool m_ShouldApproximateLoops;

      //! whether to determine the norm factor with regula falsi (true) or pythagorean approximation (false)
      bool m_DetermineAnalyticNormFactor;

      //! Excluded Volume Parameter - user provided input to adjust excluded volume in the form factor calculation
      float m_ExcludedVolumeParameter;

      //! Hydration Shell Parameter - user provided input to adjust the hydration shell in the form factor calculation
      float m_HydrationShellParameter;

      //! bool value to control side chain approximation functionality
      bool m_ShouldApproximateSideChains;

      //! Shared pointer to proposed reduced experimental data set
      util::ShPtr< storage::Vector< restraint::SasScatteringPoint> > m_ReducedExpData;

      //! euclidean distance object
      EuclideanDistance< float> m_Distance;

      //! command queue
      CommandQueue m_Queue;

      //! cl program
      cl::Program m_Program;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SaxsDebye();

      //! @brief Constructor from members
      //! @param LOOPS bool value to represent loops that are not present in the protein model
      //! @param USE_REGULA_FALSI_APPROXIMATION whether to determine the norm factor with regula falsi (true) or
      //!        pythagorean approximation (false)
      //! @param EXCLUDED_VOLUME_PARAMETER value to tune excluded volume in form factor calculation
      //! @param HYDRATION_SHELL_PARAMETER value to tune hydration shell in form factor calculation
      //! @param REDUCED_EXP_DATA shared pointer to subset of experimental data file
      SaxsDebye
      (
        const CommandQueue &QUEUE,
        const bool LOOPS = false,
        const bool USE_REGULA_FALSI_APPROXIMATION = false,
        float EXCLUDED_VOLUME_PARAMETER = 1.0,
        float HYDRATION_SHELL_PARAMETER = 0.0,
        const bool SIDE_CHAIN_APPROXIMATION = true,
        util::ShPtr< storage::Vector< restraint::SasScatteringPoint> > REDUCED_EXP_DATA =
        util::ShPtr< storage::Vector< restraint::SasScatteringPoint> >()
      );

      //! @brief Clone function
      //! @return pointer to new SaxsDebye
      SaxsDebye *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief overloaded () operator to calculate Intensity from  Q
      //! @return returns Intensity for given Q
      restraint::SasExperimentalAndCalculatedData operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write errors out to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    private:

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        UpdateQueue( GetTools());
        return io::ValidationResult( true);
      }

      //! @brief responsible for updating to a valid queue
      //! @param TOOLS opencl tools
      void UpdateQueue( Tools &TOOLS);

    }; // class SaxsDebye

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_SAXS_DEBYE_H_
