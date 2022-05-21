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
#include "nmr/bcl_nmr_rdc_container.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RDCContainer::s_Instance
    (
      GetObjectInstances().AddInstance( new RDCContainer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDCContainer::RDCContainer() :
      m_ExperimentalValues(),
      m_CalculatedValues()
    {
    }

    //! @brief construct from RDC values
    //! @param EXP_VALUES experimental values
    //! @param CALC_VALUES calculated values
    RDCContainer::RDCContainer
    (
      const storage::Vector< double> &EXP_VALUES,
      const storage::Vector< double> &CALC_VALUES
    ) :
      m_ExperimentalValues( EXP_VALUES),
      m_CalculatedValues( CALC_VALUES)
    {
      BCL_Assert
      (
        m_ExperimentalValues.GetSize() == m_CalculatedValues.GetSize(),
        "# of experimental RDS (" + util::Format()( m_ExperimentalValues.GetSize()) + ") != # of calculated RDCs(" +
          util::Format()( m_CalculatedValues.GetSize()) + ")"
      );
    }

    //! @brief Clone function
    //! @return pointer to new RDCContainer
    RDCContainer *RDCContainer::Clone() const
    {
      return new RDCContainer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDCContainer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDCContainer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExperimentalValues, ISTREAM);
      io::Serialize::Read( m_CalculatedValues, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDCContainer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExperimentalValues, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CalculatedValues, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace nmr
  
} // namespace bcl
