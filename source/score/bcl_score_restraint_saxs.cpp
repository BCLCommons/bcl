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
#include "restraint/bcl_restraint_sas_transformation.h"
#include "score/bcl_score_restraint_saxs.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RestraintSaxs::s_Instance
    (
      GetObjectInstances().AddInstance( new RestraintSaxs())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintSaxs::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "saxs_restraint");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintSaxs::RestraintSaxs() :
      m_Calc(),
      m_Score(),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief construct from member variables
    //! @param DATA_FROM_STRUCTURE_CALCULATOR function interface for calculating theoretical Saxs Curves
    //! @param DATA_AGREEMENT_CALCULATOR function interface for calculating a Saxs Curve agreement value
    //! @param SCHEME Scheme to be used
    RestraintSaxs::RestraintSaxs
    (
      const util::Implementation< restraint::SasDebyeInterface> &DEBYE_IMPLEMENTATION,
      const SasType &SCORE,
      const std::string &SCHEME
    ) :
      m_Calc( DEBYE_IMPLEMENTATION),
      m_Score( SCORE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RestraintSaxs
    RestraintSaxs *RestraintSaxs::Clone() const
    {
      return new RestraintSaxs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintSaxs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RestraintSaxs::GetAlias() const
    {
      static const std::string s_Name( "ScoreSaxs");
      return s_Name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintSaxs::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        " Computes the score between two SAXS Profiles"
      );
      parameters.AddInitializer
      (
        "",
        "algorithm to calculate the saxs profile; opencl version will run on GPU, if available",
        io::Serialization::GetAgent( &m_Calc)
      );
      parameters.Merge( m_Score.GetSerializer());
      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator() which takes a ProteinModel for calculating its agreement with the Saxs Data
    //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the Saxs Data
    //! @return return a double which is the score of the agreement of the ProteinModel with the Saxs data
    double RestraintSaxs::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // Compute the raw saxs profile
      restraint::SasExperimentalAndCalculatedData raw_data( m_Calc->operator()( PROTEIN_MODEL));

      // transform the profile into the desired form
      restraint::SasExperimentalAndCalculatedData transformed_data( restraint::SasTransformation()( raw_data));

      // calculate the score
      return m_Score( transformed_data);
    }

  } // namespace score

} // namespace bcl
