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

#ifndef BCL_RESTRAINT_HANDLER_ATOM_DISTANCE_ASSIGNED_H_
#define BCL_RESTRAINT_HANDLER_ATOM_DISTANCE_ASSIGNED_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_handler_base.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerAtomDistanceAssigned
    //! @brief HandlerAtomDistanceAssigned is for creating AtomDistance restraints from a file
    //! @details no
    //!
    //! @see @link example_restraint_handler_atom_distance_assigned.cpp @endlink
    //! @author alexanns, mendenjl
    //! @date 05/12/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API HandlerAtomDistanceAssigned :
      public HandlerBase< util::ShPtrVector< AtomDistance> >
    {

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      std::string m_Format;        //!< Format string; can be given by user
      std::string m_DefaultFormat; //!< Default format; default format for this file type
      double      m_DefaultLowerBound; //!< Default value used for lower bound if it is not contact-specific
      double      m_DefaultUpperBound; //!< Default value used for upper bound if it is not contact-specific
      double      m_DefaultDistance;   //!< Default value used for distance if it is not contact-specific

    public:

      //! @brief gives the identifying string at the top of the restraint file
      //! @return string that identifies the file as an atom distance restraint file
      static const std::string GetFileHeader();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from defaults
      HandlerAtomDistanceAssigned
      (
        const std::string &DEFAULT_EXTENSION = "",
        const std::string &DEFAULT_FORMAT = "",
        const double &LOWER_BOUND = 0.0,
        const double &UPPER_BOUND = 10.0,
        const double &DISTANCE = 8.0
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new HandlerAtomDistanceAssigned
      HandlerAtomDistanceAssigned *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads atom distance restraints from an input stream
      //! @brief input stream to read the restraints from
      //! @return the read in restraints
      util::ShPtrVector< AtomDistance> ReadRestraints( std::istream &ISTREAM) const;

    public:

      //! @brief WriteRestraints writes restraint information to an ostream
      //! TODO use list of atom distance object
      //! @param OSTREAM the stream to which the restraint information will be written
      //! @param RESTRAINT_LIST the list of restraint information that will be written to OSTREAM
      //!        the two triplets have the two chains, seqids, and atoms needed to specify the objects of the restraint
      //!        the three double in the vectornd<3> has the distance, upper bound, and lower bound, respectively
      //! @return ostream
      std::ostream &WriteRestraints
      (
        std::ostream &OSTREAM,
        const storage::Vector
        <
          storage::Pair
          <
            storage::VectorND< 2, storage::Triplet< char, int, biol::Atom> >, storage::VectorND< 3, double>
          >
        > &RESTRAINT_LIST
      ) const;

      //! @brief WriteRestraints writes restraint information to an ostream
      //! @param OSTREAM the stream to which the restraint information will be written
      //! @param RESTRAINT_LIST the list of restraints information that will be written to OSTREAM
      //! @param INCLUDE_ATOM_TYPE true to write the atom types
      //! @param INCLUDE_AA_TYPE true to write aa types
      //! @return ostream
      static std::ostream &WriteRestraints
      (
        std::ostream &OSTREAM,
        const util::ShPtrVector< AtomDistance> &RESTRAINT_LIST,
        const bool &INCLUDE_ATOM_TYPE = true,
        const bool &INCLUDE_AA_TYPE = false
      );

      //! @brief creates atom distance restraints given a set of data with distances calculated from a model ensemble
      //! @param ENSEMBLE the protein model ensemble from which distances for the data pairs will be calculated
      //! @param DATA_PAIRS the list of restraints whose distances will be calculated from the model
      //! @return list of atom distance restraints that were calculated from the data pairs and the model ensemble
      static util::ShPtrVector< AtomDistance> CreateRestraints
      (
        const assemble::ProteinEnsemble &ENSEMBLE,
        const DataSetPairwise &DATA_PAIRS
      );

      //! @brief writes a list of restraint information in rosetta format
      //! @param RESTRAINTS the list of restraint information
      //! @param OSTREAM the stream to which the restraint will be written
      //! @return std::ostream which was passed as parameter
      static std::ostream &WriteDistanceRestraintsRosettaFormat
      (
        std::ostream &OSTREAM,
        const util::ShPtrVector< AtomDistance> &RESTRAINT_LIST
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; //class HandlerAtomDistanceAssigned

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_HANDLER_ATOM_DISTANCE_ASSIGNED_H_
