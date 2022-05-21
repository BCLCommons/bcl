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

#ifndef BCL_NMR_STAR_NOE_HANDLER_H_
#define BCL_NMR_STAR_NOE_HANDLER_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_handler_base.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StarNOEHandler
    //! @brief Handler for reading/writing NOE restraints in Star format
    //! @brief Reads in NMR-Star formatted files containing NMR restraint data.  These files are
    //!        automatically generated for NMR structures deposited to the PDB and are stored at the NMR restraints grid
    //!        operated by the BMRB.  More information regarding the file format can be found at
    //!        http://www.bmrb.wisc.edu/dictionary/3.1html/SuperGroupPage.html
    //!
    //! @see @link example_nmr_star_noe_handler.cpp @endlink
    //! @author weinerbe
    //! @date Aug 9, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StarNOEHandler :
      public restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> >
    {

    private:

    //////////
    // data //
    //////////

      //! defines the lower limit for the NOEs
      static const double s_LowerLimit;

      //! defines the amount to be added to the NOE restraint for defining the upper limit
      static const double s_UpperLimit;

      //! determines the closest distance within the sequence two residues can be and still calculate a restraint
      size_t m_SequenceDistance;

      //! determines if the restraints seq ids need to be recalculated
      size_t m_SequenceOffset;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from the sequence distance and offset
      //! @param SEQUENCE_DISTANCE size_t which is the smallest distance in sequence two residues can be if a restraint
      //!        is going to be stored
      //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
      StarNOEHandler( const std::string &DEFAULT_EXTENSION = ".star", const size_t SEQUENCE_DISTANCE = 0, const size_t SEQUENCE_OFFSET = 0);

      //! @brief Clone function
      //! @return pointer to new StarNOEHandler
      StarNOEHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns short name for the class used over the command line
      //! @return short name for the class used over the command line
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads restraints from a stream
      //! @param ISTREAM the stream the restraints will be read from
      //! @return stream that the restraints were read from
      util::ShPtrVector< restraint::AtomDistance> ReadRestraints( std::istream &ISTREAM) const;

      //! @brief writes restraints to a stream
      //! @param OSTREAM the stream the restraints will be written to
      //! @param RESTRAINT the restraint that will be written to the stream
      //! @return stream the restraints were written to
      static std::ostream &WriteRestraints
      (
        std::ostream &OSTREAM,
        const util::ShPtrVector< restraint::AtomDistance> &RESTRAINT
      );

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief gets atom pair data used for determining distance restraints
      //! @param LINE_DATA the data needed to get NOE values
      //! @return vector of atom pair information (chain id, seq_id, atom type string)
      storage::Map
      <
        size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> >
      > GetAtomPairData
      (
        storage::Map
        <
          StarTagCategory,
          storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
        > &LINE_DATA
      ) const;

    }; // class StarNOEHandler

  } // namespace nmr
} // namespace bcl

#endif // BCL_NMR_STAR_NOE_HANDLER_H_
