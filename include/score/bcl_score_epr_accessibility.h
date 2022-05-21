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

#ifndef BCL_SCORE_EPR_ACCESSIBILITY_H_
#define BCL_SCORE_EPR_ACCESSIBILITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_accessibility_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EPRAccessibility
    //! @brief TODO: add a brief comment
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_score_epr_accessibility.cpp @endlink
    //! @author alexanns, rouvelgh, woetzen, karakam
    //! @date July 13, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EPRAccessibility :
      public math::FunctionInterfaceSerializable
      <
        restraint::Assignment
        <
          storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>, double, biol::AABase
        >,
        double
      >
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! map of energy functions to be used
      math::CubicSplineDamped m_EnergyFunction;

      //! the lower bounds of the bins of the histogram
      double m_HistogramLowerBound;

      //! the upper bound of the bins of the histogram
      double m_HistogramUpperBound;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are not set up from the labels
      //! @return parameters for member data that are not set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param SCHEME scheme to be used
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      EPRAccessibility
      (
        const std::string &SCHEME = GetDefaultScheme(),
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename()
      );

      //! @brief Clone function
      //! @return pointer to new EPRAccessibility
      EPRAccessibility *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns filename of the histogram being used
      //! @return filename of the histogram being used
      const std::string &GetHistogramFilename() const
      {
        return m_HistogramFileName;
      }

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief access to the energy function
      //! @return cubic spline which is "m_EnergyFunction"
      const math::CubicSplineDamped &GetEnergyFunction() const
      {
        return m_EnergyFunction;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that scores assignment of interest using pre-generated statistics
      //! @param ASSIGNMENT particular assignment of interest
      //! @return score
      double operator()
      (
        const restraint::Assignment
        <
          storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>, double, biol::AABase
        > &ASSIGNMENT
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param ASSIGNMENT
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const restraint::Assignment
        <
          storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double>, double, biol::AABase
        > &ASSIGNMENT,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @return ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    private:

      //! @brief read energy distribution for scoring accessility of amino acids based on EPR accessibility
      void ReadEnergyVector();

    }; // class EPRAccessibility

  } // namespace score

} // namespace bcl

#endif // BCL_SCORE_EPR_ACCESSIBILITY_H_
