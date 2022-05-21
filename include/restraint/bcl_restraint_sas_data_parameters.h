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

#ifndef BCL_RESTRAINT_SAS_DATA_PARAMETERS_H_
#define BCL_RESTRAINT_SAS_DATA_PARAMETERS_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasDataParameters
    //! @brief Storage Class for Saxs Parameters
    //! @details Parameters are necessary for computing form factors with a water layer and excluded volume coefficient
    //!
    //! @see @link example_restraint_saxs_data_parameters.cpp @endlink
    //! @author putnamdk
    //! @date Oct 17, 2014
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasDataParameters :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Use Sans, if True use SANS, otherwise use SAXS
      bool m_UseSansImplementation;

      //! momentum transfer point
      double m_Qvalue;

      //! Solvent Accessible Surface Area
      double m_Sasa;

      //! C1 tuning parameter for Saxs Fitting
      double m_ExcludedVolume;

      //! C2 Water layer tuning parameter for Saxs Fitting
      double m_HydrationShell;

      //! Parameter 0 < Y < 1 to represent the Deuteration level of the solvent
      double m_DeuteriumExchangeRate;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasDataParameters();

      //! @brief constructor from Q_Value alone
      explicit SasDataParameters( const double &Q_VALUE);

      //! @brief constructor from given input data
      explicit SasDataParameters
      (
        const bool &USE_SANS,
        const double &Q_VALUE,
        const double &SASA_VALUE,
        const double &EXCLUDED_VOLUME,
        const double &HYDRATION_SHELL,
        const double &DEUTERIUM_EXCHANGE_RATE
      );

      //! @brief Clone function
      //! @return pointer to new SasScatteringData
      SasDataParameters *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::str
      const std::string &GetClassIdentifier() const;

      const double &GetQValue() const
      {
        return m_Qvalue;
      }

      const double &GetSasaValue() const
      {
        return m_Sasa;
      }

      const double &GetExcludedVolume() const
      {
        return m_ExcludedVolume;
      }

      const double &GetHydrationShell() const
      {
        return m_HydrationShell;
      }

      const bool &GetSansImplementation() const
      {
        return m_UseSansImplementation;
      }

      const double &GetDeuteriumExchangeRate() const
      {
        return m_DeuteriumExchangeRate;
      }

    /////////////////
    // operations  //
    /////////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SasScatteringData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_DATA_PARAMETERS_H_
