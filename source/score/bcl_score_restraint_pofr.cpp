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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"

// includes from bcl - sorted alphabetically
#include "score/bcl_score_restraint_pofr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RestraintPofr::s_Instance
    (
      GetObjectInstances().AddInstance( new RestraintPofr())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintPofr::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "pofr_restraint");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintPofr::RestraintPofr() :
      m_Calc(),
      m_Score(),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief construct from member variables
    //! @param DATA_FROM_STRUCTURE_CALCULATOR function interface for calculating theoretical Saxs Curves
    //! @param DATA_AGREEMENT_CALCULATOR function interface for calculating a Saxs Curve agreement value
    //! @param SCHEME Scheme to be used
    RestraintPofr::RestraintPofr
    (
      const util::Implementation< restraint::SasPofRInterface> &POFR_IMPLEMENTATION,
      const PofR &SCORE,
      const std::string &SCHEME
    ) :
      m_Calc( POFR_IMPLEMENTATION),
      m_Score( SCORE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RestraintPofr
    RestraintPofr *RestraintPofr::Clone() const
    {
      return new RestraintPofr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintPofr::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RestraintPofr::GetAlias() const
    {
      static const std::string s_Name( "ScorePofr");
      return s_Name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintPofr::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        " Computes the difference between experimental SAS and computed density histograms"
      );
      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator() which takes a ProteinModel for calculating its agreement with the Pofr Data
    //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the Pofr Data
    //! @return return a double which is the score of the agreement of the ProteinModel with the Pofr data
    double RestraintPofr::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      //BCL_Message( util::Message::e_Standard, "Inside Operator to compute profile and score agreement");

      // Compute the raw saxs profile
      // m_Calc contains the experimental data stored as a data member
      restraint::SasExperimentalAndCalculatedDensity raw_data( m_Calc->operator()( PROTEIN_MODEL));

      // calculate the score
      return m_Score( raw_data);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RestraintPofr::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RestraintPofr::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace score

} // namespace bcl
