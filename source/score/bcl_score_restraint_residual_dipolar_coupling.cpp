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
#include "score/bcl_score_restraint_residual_dipolar_coupling.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "nmr/bcl_nmr_rdc_container.h"
#include "restraint/bcl_restraint_rdc.h"
#include "restraint/bcl_restraint_rdc_assignment.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> RestraintResidualDipolarCoupling::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new RestraintResidualDipolarCoupling())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintResidualDipolarCoupling::RestraintResidualDipolarCoupling() :
      m_CalculatorOfRDCs(),
      m_RDCAgreementCalculator(),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief construct from member variables
    //! @param RDCS experimental rdcs
    //! @param RDC_FROM_STRUCTURE_CALCULATOR function interface for calculating theoretical RDCs
    //! @param RDC_AGREEMENT_CALCULATOR function interface for calculating an RDC agreement value
    //! @param SCHEME Scheme to be used
    RestraintResidualDipolarCoupling::RestraintResidualDipolarCoupling
    (
      const util::ShPtr< restraint::RDC> &RDCS,
      const math::FunctionInterfaceSerializable< restraint::RDCAssignment, nmr::RDCContainer> &RDC_FROM_STRUCTURE_CALCULATOR,
      const math::FunctionInterfaceSerializable< nmr::RDCContainer, double> &RDC_AGREEMENT_CALCULATOR,
      const std::string &SCHEME
    ) :
      m_CalculatorOfRDCs( RDC_FROM_STRUCTURE_CALCULATOR),
      m_RDCAgreementCalculator( RDC_AGREEMENT_CALCULATOR),
      m_RDCs( RDCS),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RestraintResidualDipolarCoupling
    RestraintResidualDipolarCoupling *RestraintResidualDipolarCoupling::Clone() const
    {
      return new RestraintResidualDipolarCoupling( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RestraintResidualDipolarCoupling::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintResidualDipolarCoupling::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "rdc_restraint");

      // end
      return s_default_scheme;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &RestraintResidualDipolarCoupling::GetAlias() const
    {
      static const std::string s_name( "RestraintResidualDipolarCoupling");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintResidualDipolarCoupling::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Scores a protein model's agreement with residual dipolar coupling data"
      );
      serializer.AddInitializer
      (
        "calculator of RDCs",
        "Takes a list of assigned RDCs and calculates theoretical RDCs",
        io::Serialization::GetAgent( &m_CalculatorOfRDCs)
      );
      serializer.AddInitializer
      (
        "RDC agreement calculator",
        "takes a ResidualDipolarCouplingContainer and returns an agreement value",
        io::Serialization::GetAgent( &m_RDCAgreementCalculator)
      );
      return serializer;
    }
  ////////////////
  // operations //
  ////////////////

    //! @brief operator() which takes an ProteinModel for calculating its agreement with the RDC data
    //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the RDC data
    //! @return return a double which is the score of the agreement of the ProteinModel with the RDC data
    double RestraintResidualDipolarCoupling::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // if no restraints were found
      if( !m_RDCs.IsDefined())
      {
        // return a score of zero
        BCL_MessageStd( "No RDC restraints found, not scoring the model");
        return 0.0;
      }

      // calculate the RDCs and the score
      const double score
      (
        m_RDCAgreementCalculator->operator ()
        (
          m_CalculatorOfRDCs->operator ()( m_RDCs->GenerateAssignment( PROTEIN_MODEL))
        )
      );

      // weight score by log(#restraints + 1) * (#AAs)
      return score *
          log10( double( m_RDCs->GetData().GetSize() + 1)) *
          double( PROTEIN_MODEL.GetNumberAAs());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RestraintResidualDipolarCoupling::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_RDCs, ISTREAM);
      io::Serialize::Read( m_CalculatorOfRDCs, ISTREAM);
      io::Serialize::Read( m_RDCAgreementCalculator, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RestraintResidualDipolarCoupling::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_RDCs, OSTREAM, INDENT);
      io::Serialize::Write( m_CalculatorOfRDCs, OSTREAM, INDENT);
      io::Serialize::Write( m_RDCAgreementCalculator, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the RDC data
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &RestraintResidualDipolarCoupling::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      // if no restraints were found
      if( !m_RDCs.IsDefined())
      {
        return OSTREAM;
      }

      // call the WriteDetailedSchemeAndValues from the agreement calculator
      return m_RDCAgreementCalculator->WriteDetailedSchemeAndValues
          (
            m_CalculatorOfRDCs->operator ()( m_RDCs->GenerateAssignment( PROTEIN_MODEL)),
            OSTREAM
          );
    }

  } // namespace score
} // namespace bcl
