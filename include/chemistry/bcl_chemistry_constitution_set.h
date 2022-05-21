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

#ifndef BCL_CHEMISTRY_CONSTITUTION_SET_H_
#define BCL_CHEMISTRY_CONSTITUTION_SET_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_hydrogens_handler.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstitutionSet
    //! @brief Container class for FragmentConstitutionShared objects
    //!
    //! @see @link example_chemistry_constitution_set.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 25, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConstitutionSet :
      public util::ObjectInterface
    {
    public:

    //////////////
    // typedefs //
    //////////////

      typedef util::ShPtrList< FragmentConstitutionShared>::const_iterator const_iterator;

    private:

    //////////
    // data //
    //////////

      //! list of FragmentConstitutionShared
      util::ShPtrList< FragmentConstitutionShared> m_Constitutions;

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Node
      //! @brief Storage class for information about a molecule in the set
      //!
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Oct 24, 2012
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Node
      {
        const_iterator                     m_Itr;   //!< iterator on m_Configurations
        graph::ConstGraph< size_t, size_t> m_Graph; //!< ConstGraph of the molecule

        //! @brief test whether the graph is cached
        //! @return true if the graph is cached on this node
        bool IsGraphCached() const
        {
          return m_Graph.GetSize() > size_t( 0);
        }

        //! @brief cache the graph for the molecule using the iterator
        void CacheGraph();

      };

      //! map of hash string to itr and graph::ConstGraph
      std::map< std::string, std::list< Node> > m_HashToNodeMap;

      //! an isomorphism comparer, cached to avoid constantly reallocating memory
      //! mutable because it will be changed in Find, even though the set itself does not logically change
      mutable graph::SubgraphIsomorphism< size_t, size_t> m_Isomorphism;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConstitutionSet()
      {}

      //! @brief construct from iterator of ConstitutionInterface
      //! @param MOLECULES an iterator to configuration of molecules
      ConstitutionSet( iterate::Generic< const ConstitutionInterface> MOLECULES);

      //! @brief Clone function
      //! @return pointer to new ConstitutionSet
      ConstitutionSet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return a list containing the constitutions
      //! @return a list containing the constitutions
      const util::ShPtrList< FragmentConstitutionShared> &GetConstitutions() const
      {
        return m_Constitutions;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the number of constitutions
      //! @return the number of constitutions
      size_t GetSize() const
      {
        return m_Constitutions.GetSize();
      }

      //! @brief returns const iterator begin to constitutions
      //! @return const iterator begin to constitutions
      const_iterator Begin() const
      {
        return m_Constitutions.Begin();
      }

      //! @brief returns const iterator end to constitutions
      //! @return const iterator end to constitutions
      const_iterator End() const
      {
        return m_Constitutions.End();
      }

      //! @brief find a particular constitution
      //! @param FRAGMENT a constitution to search for
      //! @param ISOMORPHISM isomorphism vector to store isomorphism
      //! @return an iterator to the constitution in this set
      const_iterator Find
      (
        const ConstitutionInterface &FRAGMENT,
        util::SiPtr< storage::Vector< size_t> > ISOMORPHISM = util::SiPtr< storage::Vector< size_t> >()
      ) const;

      //! @brief append a molecule to SetConstitution
      //! @param FRAGMENT fragment constitution shared that needs to be added to set constitution
      //! @param ISOMORPHISM isomorphism vector to store isomorphism
      //! @return a pair containing iterator to constitution layer of FRAGMENT and whether constitution has been seen
      //!         earlier (false) or seen for the first time (true)
      std::pair< const_iterator, bool> Insert
      (
        const ConstitutionInterface &FRAGMENT,
        util::SiPtr< storage::Vector< size_t> > ISOMORPHISM = util::SiPtr< storage::Vector< size_t> >()
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief make a hash string of the atom types and bond types
      //! @param a fragment constitution
      //! @return the hash string of atom types and bond types for the given fragment
      static std::string MakeBasicHashString( const ConstitutionInterface &FRAGMENT);

      //! @brief function used internally to search for a given fragment within member data structures
      //! @param FRAGMENT a constitution to search for
      //! @param NODES molecules to compare FRAGMENT against
      //! @param ISOMORPHISM isomorphism vector to store isomorphism if it is calculated
      //! @param GRAPH reference to a const graph, will be set to the graph for FRAGMENT if a graph is created
      const_iterator InternalFind
      (
        const ConstitutionInterface &FRAGMENT,
        const std::list< Node> &NODES,
        storage::Vector< size_t> &ISOMORPHISM,
        graph::ConstGraph< size_t, size_t> &GRAPH
      ) const;

    }; // class ConstitutionSet

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONSTITUTION_SET_H_
