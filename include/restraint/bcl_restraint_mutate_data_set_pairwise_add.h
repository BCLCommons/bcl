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

#ifndef BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_ADD_H_
#define BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_ADD_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_data_set_pairwise.h"
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDataSetPairwiseAdd
    //! @brief takes a DataSetPairwise and mutates it by adding some number of DataPairs from an eligable pool
    //! @details Contains the pool of all possible DataPairs that can be added as well as a range of possible number
    //!          of DataPairs to add during a singe mutation.
    //!
    //! @see @link example_restraint_mutate_data_set_pairwise_add.cpp @endlink
    //! @author alexanns
    //! @date May 7, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateDataSetPairwiseAdd :
      public math::MutateInterface< DataSetPairwise>
    {

    private:

    //////////
    // data //
    //////////

      //! the pool of possible DataPairwise objects that can be added
      util::ShPtr< DataSetPairwise> m_CompleteDataSet;

      //! the range of number of DataPairwise objects that can be added
      math::Range< size_t> m_SizeRange;

      //! scheme of this mutate
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from optional scheme
      //! @param SCHEME optional scheme
      explicit MutateDataSetPairwiseAdd( const std::string &SCHEME = GetDefaultScheme());

      //! @brief constructor taking sequence
      //! @param POOL_DATA_SET pool of possible data points to add
      //! @param SCHEME optional scheme
      explicit MutateDataSetPairwiseAdd
      (
        const util::ShPtr< DataSetPairwise> &POOL_DATA_SET, const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief constructor taking sequence
      //! @param POOL_DATA_SET pool of possible data points to add
      //! @param MIN min possible number that will be added
      //! @param MAX max possible number that will be added
      //! @param SCHEME optional scheme
      MutateDataSetPairwiseAdd
      (
        const util::ShPtr< DataSetPairwise> &POOL_DATA_SET, const size_t MIN, const size_t MAX,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDataSetPairwiseAdd
      MutateDataSetPairwiseAdd *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
      //! @param DATA_SET DataSetPairwise of interest that will be mutated
      //! @return MutateResult that results from mutating to the argument DATA_SET
      math::MutateResult< DataSetPairwise> operator()( const DataSetPairwise &DATA_SET) const;

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

    public:

    }; // class MutateDataSetPairwiseAdd

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_MUTATE_DATA_SET_PAIRWISE_ADD_H_ 
