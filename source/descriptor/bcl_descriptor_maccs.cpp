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
#include "descriptor/bcl_descriptor_maccs.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_ifstream.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MACCS::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MACCS()
      )
    );

    namespace
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class BondCheck
      //! @brief simple class to compare bond types; allowing type 5 (~) to match any bond, and type 4 (aromatic)
      //!        to match single, double, and aromatic bonds
      //! @remarks example unnecessary
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct BondCheck :
        public std::binary_function< size_t, size_t, bool>
      {
        //! @brief comparison for the bond types
        bool operator()( const size_t &ARGUMENT1, const size_t &ARGUMENT2) const
        {
          // allow type 5 (~) to match any bond, and type 4 (aromatic)
          // to match single, double, and aromatic bonds
          return ARGUMENT1 == ARGUMENT2
                 || ( ARGUMENT2 == 5 && ARGUMENT1 < 4)
                 || ( ARGUMENT1 == 4 && ARGUMENT2 < 5 && ARGUMENT2 != 3);
        }
      };

      //! @brief get a map from pair of element types to position in the maccs key
      storage::Map< std::pair< chemistry::ElementType, chemistry::ElementType>, int> SimpleBondElementTypesMap()
      {
        static storage::Map< std::pair< chemistry::ElementType, chemistry::ElementType>, int> s_map;
        if( s_map.IsEmpty())
        {
          const chemistry::ElementTypes &elements( chemistry::GetElementTypes());
          ///Lithium/////////////////////////////////////////////////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Hydrogen)]          = 0;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Lithium)]           = 1;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Boron)]             = 2;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Carbon)]            = 3;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Oxygen)]            = 4;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Fluorine)]          = 5;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Phosphorus)]        = 6;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Sulfur)]            = 7;
          s_map[ std::make_pair( elements.e_Lithium, elements.e_Chlorine)]          = 8;
          //////////Boron////////////////////////////////////////////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Boron, elements.e_Hydrogen)]            = 9;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Boron)]               = 10;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Carbon)]              = 11;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Nitrogen)]            = 12;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Oxygen)]              = 13;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Fluorine)]            = 14;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Silicon)]             = 15;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Phosphorus)]          = 16;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Sulfur)]              = 17;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Chlorine)]            = 18;
          s_map[ std::make_pair( elements.e_Boron, elements.e_Bromine)]             = 19;
          /////////////////Carbon////////////////////////////////////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Hydrogen)]           = 20;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Carbon)]             = 21;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Nitrogen)]           = 22;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Oxygen)]             = 23;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Fluorine)]           = 24;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Sodium)]             = 25;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Magnesium)]          = 26;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Aluminum)]           = 27;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Silicon)]            = 28;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Phosphorus)]         = 29;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Sulfur)]             = 30;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Chlorine)]           = 31;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Arsenic)]            = 32;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Selenium)]           = 33;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Bromine)]            = 34;
          s_map[ std::make_pair( elements.e_Carbon, elements.e_Iodine)]             = 35;
          ///////////////////////Nitrogen////////////////////////////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Hydrogen)]         = 36;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Nitrogen)]         = 37;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Oxygen)]           = 38;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Fluorine)]         = 39;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Silicon)]          = 40;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Phosphorus)]       = 41;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Sulfur)]           = 42;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Chlorine)]         = 43;
          s_map[ std::make_pair( elements.e_Nitrogen, elements.e_Bromine)]          = 44;
          ///////////////////////////////Oxygen//////////////////////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Hydrogen)]           = 45;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Oxygen)]             = 46;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Magnesium)]          = 47;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Sodium)]             = 48;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Aluminum)]           = 49;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Silicon)]            = 50;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Phosphorus)]         = 51;
          s_map[ std::make_pair( elements.e_Oxygen, elements.e_Potassium)]          = 52;
          /////////////////////////////////////Fluorine//////////////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Fluorine, elements.e_Phosphorus)]       = 53;
          s_map[ std::make_pair( elements.e_Fluorine, elements.e_Sulfur)]           = 54;
          /////////////////////////////////////////////Silicon///////////////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Silicon, elements.e_Hydrogen)]          = 55;
          s_map[ std::make_pair( elements.e_Silicon, elements.e_Silicon)]           = 56;
          s_map[ std::make_pair( elements.e_Silicon, elements.e_Chlorine)]          = 57;
          /////////////////////////////////////////////////////Arsenic///////////////////////////////////////////////////
          s_map[ std::make_pair( elements.e_Arsenic, elements.e_Hydrogen)]          = 58;
          s_map[ std::make_pair( elements.e_Arsenic, elements.e_Arsenic)]           = 59;
        }

        return s_map;
      }

      //! @brief create graphs for all the smarts strings at the end of this file
      storage::Vector< graph::ConstGraph< size_t, size_t> > SmartsStringsInputGraphs
      (
        const std::string SMARTS[]
      )
      {
        storage::Vector< graph::ConstGraph< size_t, size_t> > graphs;

        const size_t undef_bond_type( util::GetUndefined< size_t>());
        const std::string bond_type_string( "!-=#:~");
        for( size_t i( 0); i < size_t( MACCS::s_NumberSMARTSQueries); ++i)
        {
          storage::Vector< size_t> connection_stack;
          storage::Vector< size_t> element_vector;
          storage::Vector< size_t> aromatic_stack;
          storage::Vector< graph::UndirectedEdge< size_t> > bond_vector;
          size_t bond_type( undef_bond_type);
          size_t ring_bond_id( util::GetUndefined< size_t>());

          for
          (
            std::string::const_iterator
              string_itr( SMARTS[ i].begin()), string_itr_end( SMARTS[ i].end());
            string_itr != string_itr_end;
            ++string_itr
          )
          {
            if( std::isalpha( *string_itr))
            {
              bool aromatic( false);
              std::string element_symbol( 1, *string_itr);
              if( std::islower( *string_itr))
              {
                aromatic = true;
                element_symbol[ 0] = std::toupper( element_symbol[ 0]);
              }

              // checks for second letter example: Cl
              if( !aromatic && ( string_itr + 1) != string_itr_end && std::islower( *( string_itr + 1)))
              {
                element_symbol += *++string_itr;    //skips the "l" char in the next iteration
              }

              // get the atomic number
              size_t element_type_lookup( chemistry::GetElementTypes().ElementTypeLookup( element_symbol).GetIndex());
              if( !element_vector.IsEmpty()) // is this any atom other than the first
              {
                // yes so assign bond
                bond_vector.PushBack
                (
                  graph::UndirectedEdge< size_t>
                  (
                    connection_stack.LastElement(),
                    element_vector.GetSize(),
                    util::IsDefined( bond_type) ? bond_type : ( aromatic ? 4 : 1)
                  )
                );
                connection_stack.LastElement() = element_vector.GetSize();
                aromatic_stack.LastElement() = aromatic;
              }
              else
              {
                //no so pushback connection stack
                connection_stack.PushBack( element_vector.GetSize());
                aromatic_stack.PushBack( aromatic);
              }
              element_vector.PushBack( element_type_lookup);
              bond_type = undef_bond_type;
            }
            else if( *string_itr == '(')    // start of a branch
            {
              // copy last element of connection_stack to new element
              connection_stack.PushBack( connection_stack.LastElement());
              aromatic_stack.PushBack( 0);
            }
            else if( *string_itr == ')')    // end of branch
            {
              //remove last element from connection stack
              connection_stack.PopBack();
              aromatic_stack.PopBack();
            }
            else if( *string_itr == '1')
            {
              if( util::IsDefined( ring_bond_id))
              {
                bond_vector.PushBack
                (
                  graph::UndirectedEdge< size_t>
                  (
                    connection_stack.LastElement(),
                    ring_bond_id,
                    util::IsDefined( bond_type) ? bond_type : ( aromatic_stack.LastElement() ? 4 : 1)
                  )
                );
                ring_bond_id = undef_bond_type;
              }
              else
              {
                ring_bond_id = element_vector.GetSize() - 1;
              }
            }
            else //Char must be a bond
            {
              bond_type = bond_type_string.find( *string_itr);
            }
          }
          // Make graph with only vertices with element vector
          graphs.PushBack( graph::ConstGraph< size_t, size_t>( element_vector, bond_vector, undef_bond_type));
        }

        return graphs;
      } //end function call

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    // rest all histograms!
    MACCS::MACCS() :
      m_ElementHistogram( chemistry::GetElementTypes().GetEnumCount(), size_t( 0)),
      m_RingSizeHistogram( 8, 0),
      m_RingTypeHistogram( 9, storage::Vector< size_t>( 8, 0))
    {
    }

    //! @brief virtual copy constructor
    MACCS *MACCS::Clone() const
    {
      return new MACCS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MACCS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MACCS::GetAlias() const
    {
      static const std::string s_name( "MACCS");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MACCS::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      //reset all histograms to zero once start for each molecule
      m_ElementHistogram.SetAllElements( 0);
      m_RingSizeHistogram.SetAllElements( 0);
      for( size_t i( 0); i < 9; ++i)
      {
        m_RingTypeHistogram( i).SetAllElements( 0);
      }

    ////////////////////////////////////////////////
    // Section 1: Check for Element Types Present //
    ////////////////////////////////////////////////

      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);    //allows data type (use) of molecule

      chemistry::ConformationGraphConverter graph_converter                 // creates a object called graph_converter
      (
        chemistry::ConformationGraphConverter::e_ElementType,
        chemistry::ConfigurationalBondTypeData::e_ConfigurationalBondType
      );

      graph::ConstGraph< size_t, size_t> graph_element_bondtype( graph_converter( molecule));

      // iterate through all atoms of the molecule and count their appearance
      for( iterate::Generic< const chemistry::AtomConformationalInterface> itr( molecule.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        if( itr->GetElementType().IsDefined())
        {
          ++m_ElementHistogram( itr->GetElementType().GetIndex());   // we have a graph and a histogram
          m_ElementHistogram( 0) += itr->GetNumberofValenceBondsWithOrder( 1);
        }
      }

      size_t bit_position( 0);
      for( size_t element_array_index( 0); element_array_index < 91; ++element_array_index)   //loop through all chemical elements
      {
        size_t initial_check( 1);
        if( element_array_index == 0)   // Hydrogen tests count at 4,8,16,32 instead of the others at 1,2,4,8 etc.
        {                               // Vector position is atomic number - 1
          initial_check = 4;
        }
        else if( element_array_index == 5)   //Carbon tests count at 2,4,8,16,32
        {
          initial_check = 2;
        }

        for( size_t current_bit( 0); current_bit < s_ElementArray[ element_array_index]; ++current_bit)
        {
          // testing to see if chemical element count is greater than 2^x
          if( m_ElementHistogram( element_array_index) >= ( initial_check << current_bit)) //access first element of vector (Vector[0])
          {
            STORAGE( bit_position + current_bit) = 1;
          }
          else
          {
            STORAGE( bit_position + current_bit) = 0;
          }
        }
        bit_position += s_ElementArray[element_array_index];    //keep sum of how many array elements used to determine placement of new array elements
      }

    ////////////////////////////////////
    // Section 2: Check for Ring sets //
    ////////////////////////////////////

      graph::EdgeCoverRingPerception edge_cover_ring( graph_element_bondtype); //constructor call for the RingPerception Class

      storage::List< graph::Ring> list_rings( edge_cover_ring.GetRings());   //returns a list <Ring> (a list of the class Ring)

      size_t carbon_index( chemistry::GetElementTypes().e_Carbon.GetIndex());
      size_t nitrogen_index( chemistry::GetElementTypes().e_Nitrogen.GetIndex());

      // determine how many carbon aromatic and hetero aromatic rings there are, all in one loop
      size_t number_carbon_aromatic_rings( 0);
      size_t number_hetero_aromatic_rings( 0);
      // iterate through list of Rings
      for
      (
          storage::List< graph::Ring>::const_iterator itr_ring( list_rings.Begin()), itr_ring_end( list_rings.End());
          itr_ring != itr_ring_end;
          ++itr_ring
      )
      {
        size_t number_aromatic_bonds( 0);
        size_t number_single_bonds( 0);
        size_t number_double_bonds( 0);
        size_t number_triple_bonds( 0);
        size_t number_nitrogen_atoms( 0);
        size_t number_carbon_atoms( 0);
        size_t ring_size( itr_ring->GetSize());

        //Iterate through Rings
        for
        (
          graph::Ring::const_iterator itr( itr_ring->Begin()), itr_prev( --itr_ring->End()), itr_end( itr_ring->End());
          itr != itr_end;
          itr_prev = itr, ++itr
        )
        {
          // get the bond type index between the atoms indicated by itr and itr_prev
          const size_t bond_type_index( graph_element_bondtype.GetEdgeData( *itr, *itr_prev));
          chemistry::ConfigurationalBondType bond_type( bond_type_index);

          const size_t bond_order_aromatic
          (
            bond_type->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic)
          );
          if( bond_order_aromatic == 1)
          {
            ++number_single_bonds;
          }
          else if( bond_order_aromatic == 2)
          {
            ++number_double_bonds;
          }
          else if( bond_order_aromatic == 3)
          {
            ++number_triple_bonds;
          }
          else if( bond_order_aromatic == 4)
          {
            ++number_aromatic_bonds;
          }
        }
        // Iterate through rings for atom type count
        for
        (
          graph::Ring::const_iterator itr_element( itr_ring->Begin()), itr_element_end( itr_ring->End());
          itr_element != itr_element_end;
          ++itr_element
        )
        {
          if( graph_element_bondtype.GetVertexData( *itr_element) == carbon_index)
          {
            ++number_carbon_atoms;
          }
          else if( graph_element_bondtype.GetVertexData( *itr_element) == nitrogen_index)
          {
            ++number_nitrogen_atoms;
          }
        }

        if( number_aromatic_bonds == ring_size && number_carbon_atoms == ring_size)
        {
          ++number_carbon_aromatic_rings;           // tally of aromatic rings in molecule
        }
        else if( number_aromatic_bonds == ring_size && ( number_carbon_atoms + number_nitrogen_atoms < ring_size))
        {
          ++number_hetero_aromatic_rings;         // tally of hetero aromatic rings
        }

        //tally values into the  m_RingTypeHistogram
        size_t storage_index( 0);
        if( ring_size < 11)
        {
          ++m_RingSizeHistogram(ring_size - 3);   // tally the ring size histogram vector

          //////////check for carbon only ring////////////////////
          if( number_carbon_atoms == ring_size)
          {
            if( number_single_bonds == ring_size) //check for saturated ring
            {
              ++m_RingTypeHistogram( storage_index)( ring_size - 3);
            }
            else if( number_aromatic_bonds == ring_size) // check for aromatic
            {
              ++m_RingTypeHistogram( storage_index + 2)( ring_size - 3);
            }
            else // check for unsaturated
            {
              ++m_RingTypeHistogram( storage_index + 1)( ring_size - 3);
            }
          }

          //////////check for any ring that has a nitrogen atom/////////////////
          if( number_nitrogen_atoms > 0)    //Check for nitrogen in ring
          {
            if( number_single_bonds == ring_size) //check for saturated ring
            {
              ++m_RingTypeHistogram( storage_index + 3)( ring_size - 3);
            }
            else if( number_aromatic_bonds == ring_size) // check for aromatic
            {
              ++m_RingTypeHistogram( storage_index + 5)( ring_size - 3);
            }
            else // check for unsaturated
            {
              ++m_RingTypeHistogram( storage_index + 4)( ring_size - 3);
            }
          }
          /////////check for hetero ring (carbon - nitrogen only rings are not counted again)///////////////////
          if( ( number_carbon_atoms + number_nitrogen_atoms) < ring_size) //check for hetero ring (nitrogen (only) heteroatoms are not counted twice)
          {
            if( number_single_bonds == ring_size) //check for saturated ring
            {
              ++m_RingTypeHistogram( storage_index + 6)( ring_size - 3);
            }
            else if( number_aromatic_bonds == ring_size) // check for aromatic
            {
              ++m_RingTypeHistogram( storage_index + 8)( ring_size - 3);
            }
            else // check for unsaturated
            {
              ++m_RingTypeHistogram( storage_index + 7)( ring_size - 3);
            }
          }

        }
      }

      //******** assign to "bit" location for ring sizes *******
      //                 Max count for rings of size 3, 4, 5, 6, 7, 8, 9, 10
      static size_t s_numer_times_ring_checked[] = { 2, 2, 5, 5, 2, 2, 1, 1};
      for( size_t ring_size_index = 0; ring_size_index < 8; ++ring_size_index)
      {
        for( size_t cur_count( 0); cur_count < s_numer_times_ring_checked[ ring_size_index]; ++cur_count, ++bit_position)
        {
          STORAGE( bit_position) = ( m_RingSizeHistogram( ring_size_index) > cur_count);
        }
      }

      //********  assign to "bit" location for aromatic and hetero-aromatic ring types *only* ******
      for( size_t ring_type_index( 0); ring_type_index < 4; ++ring_type_index)
      {
        STORAGE( bit_position++) = ( number_carbon_aromatic_rings > ring_type_index);
        STORAGE( bit_position++) = ( number_hetero_aromatic_rings > ring_type_index);
      }

      //********   assign to "bit" location of specific ring types ******

      for( size_t outer_vector( 0); outer_vector < 9; ++outer_vector)
      {
        for( size_t inner_vector( 0); inner_vector < 8; ++inner_vector)
        {
          for( size_t cur_count( 0); cur_count < s_numer_times_ring_checked[ inner_vector]; ++cur_count, ++bit_position)
          {
            STORAGE( bit_position) = ( m_RingTypeHistogram( outer_vector)( inner_vector) > cur_count);
          }
        }
      }

      //////////////////////////////////////////
      ////// Section 3: Atom Pairs//////////////
      //////////////////////////////////////////
      static const storage::Map< std::pair< chemistry::ElementType, chemistry::ElementType>, int>
        s_simple_bond_map( SimpleBondElementTypesMap());

      storage::Vector< sdf::BondInfo> bond_info( molecule.GetBondInfo());   //gets all bonds in the graph

      for
      (
        storage::Vector<sdf::BondInfo>::const_iterator itr_bond( bond_info.Begin()), itr_bond_end( bond_info.End());    //iterate through list of bonds
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {

        chemistry::ElementType low_atom( graph_element_bondtype.GetVertices()( itr_bond->GetAtomIndexLow()));
        chemistry::ElementType high_atom( graph_element_bondtype.GetVertices()( itr_bond->GetAtomIndexHigh()));

        storage::Map< std::pair< chemistry::ElementType, chemistry::ElementType>, int>::const_iterator map_itr
        (
          s_simple_bond_map.Find( std::make_pair( low_atom, high_atom))
        );

        if( map_itr == s_simple_bond_map.End())
        {
          //not found. check for reverse case
          map_itr = s_simple_bond_map.Find( std::make_pair( high_atom, low_atom));
        }
        if( map_itr != s_simple_bond_map.End())
        {
          STORAGE( map_itr->second + bit_position) = 1;
        }
      }
      //bit position at 320

      bit_position += 64;

      //bit position at 384

      //////////////////////////////////////////
      /////// Section 4, 5 and 6: SMARTS  //////
      //////////////////////////////////////////

      chemistry::ConformationGraphConverter smarts_graph_converter
      (
        chemistry::ConformationGraphConverter::e_ElementType,
        chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic
      );
      graph::ConstGraph< size_t, size_t> graph_bondorderaro_element( smarts_graph_converter( molecule));

      //s_SMARTS_const_graphs is a vector of all const_graphs
      static storage::Vector< graph::ConstGraph< size_t, size_t> > s_SMARTS_const_graphs
      (
        SmartsStringsInputGraphs( s_AroSMARTSArray)
      );
      util::ShPtr< util::BinaryFunctionInterface< size_t, size_t, bool> > bond_check_ptr
      (
        new util::BinaryFunctionSTLWrapper< BondCheck>()
      );

      util::OwnPtr< graph::ConstGraph< size_t, size_t> > graph_bondorderaro_element_ptr
      (
        &graph_bondorderaro_element,
        false
      );
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_sub_isomorphism
      (
        graph::CommonSubgraphIsomorphism< size_t, size_t>::e_Connected,
        graph::CommonSubgraphIsomorphism< size_t, size_t>::GetDefaultVertexComparison(),
        bond_check_ptr
      );

      common_sub_isomorphism.SetGraphA( graph_bondorderaro_element_ptr);
      for
      (
        storage::Vector< graph::ConstGraph< size_t, size_t> >::iterator
        itr_frags( s_SMARTS_const_graphs.Begin()), itr_frags_end( s_SMARTS_const_graphs.End());
        itr_frags != itr_frags_end;
        ++itr_frags, ++bit_position
      )
      {
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > sub_graph_ptr( &*itr_frags, false);
        common_sub_isomorphism.SetGraphB( sub_graph_ptr);
        size_t size_sub_graph( itr_frags->GetSize());

        //check for isomorphism
        common_sub_isomorphism.FindIsomorphism( common_sub_isomorphism.EstimateUpperBounds(), size_sub_graph);
        if( common_sub_isomorphism.GetIsomorphism().GetSize() == size_sub_graph)
        {
          STORAGE( bit_position) = 1; //bit placement 384 to 769
        }
      }
    } // calculate

    /////////////////////////////
    /// Helper Functions/////////
    /////////////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MACCS::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Retrieves the MACCS of each molecule using a modified pubchem MACCS "
        "fingerprint. See ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt for details on "
        "the pubchem values"
      );
      return parameters;
    }

    const std::string MACCS::s_AroSMARTSArray[ s_NumberSMARTSQueries] =
    {
      ///Section 4/// Starts at bit 384
      "C(~Br)(~C)", "C(~Br)(~C)(~C)", "C(~Br)(~H)", "C(~Br)(:C)", "C(~Br)(:N)", "C(~C)(~C)", "C(~C)(~C)(~C)",
      "C(~C)(~C)(~C)(~C)", "C(~C)(~C)(~C)(~H)", "C(~C)(~C)(~C)(~N)", "C(~C)(~C)(~C)(~O)", "C(~C)(~C)(~H)(~N)",
      "C(~C)(~C)(~H)(~O)", "C(~C)(~C)(~N)", "C(~C)(~C)(~O)", "C(~C)(~Cl)", "C(~C)(~Cl)(~H)", "C(~C)(~H)",
      "C(~C)(~H)(~N)", "C(~C)(~H)(~O)", "C(~C)(~H)(~O)(~O)", "C(~C)(~H)(~P)", "C(~C)(~H)(~S)", "C(~C)(~I)", "C(~C)(~N)",
      "C(~C)(~O)", "C(~C)(~S)", "C(~C)(~Si)", "C(~C)(:C)", "C(~C)(:C)(:C)", "C(~C)(:C)(:N)", "C(~C)(:N)",
      "C(~C)(:N)(:N)", "C(~Cl)(~Cl)", "C(~Cl)(~H)", "C(~Cl)(:C)", "C(~F)(~F)", "C(~F)(:C)", "C(~H)(~N)", "C(~H)(~O)",
      "C(~H)(~O)(~O)", "C(~H)(~S)", "C(~H)(~Si)", "C(~H)(:C)", "C(~H)(:C)(:C)", "C(~H)(:C)(:N)", "C(~H)(:N)",
      "C(~H)(~H)(~H)", "C(~N)(~N)", "C(~N)(:C)", "C(~N)(:C)(:C)", "C(~N)(:C)(:N)", "C(~N)(:N)", "C(~O)(~O)",
      "C(~O)(:C)", "C(~O)(:C)(:C)", "C(~S)(:C)", "C(:C)(:C)", "C(:C)(:C)(:C)", "C(:C)(:C)(:N)", "C(:C)(:N)",
      "C(:C)(:N)(:N)", "C(:N)(:N)", "N(~C)(~C)", "N(~C)(~C)(~C)", "N(~C)(~C)(~H)", "N(~C)(~H)", "N(~C)(~H)(~N)",
      "N(~C)(~O)", "N(~C)(:C)", "N(~C)(:C)(:C)", "N(~H)(~N)", "N(~H)(:C)", "N(~H)(:C)(:C)", "N(~O)(~O)", "N(~O)(:O)",
      "N(:C)(:C)", "N(:C)(:C)(:C)", "O(~C)(~C)", "O(~C)(~H)", "O(~C)(~P)", "O(~H)(~S)", "O(:C)(:C)", "P(~C)(~C)",
      "P(~O)(~O)", "S(~C)(~C)", "S(~C)(~H)", "S(~C)(~O)", "Si(~C)(~C)",
      ///Section 5/// // starts at bit 473
      "C=C", "C#C", "C=N", "C#N", "C=O", "C=S", "N=N", "N=O", "N=P", "P=O", "P=P",
      "C(#C)(-C)", "C(#C)(-H)", "C(#N)(-C)", "C(-C)(-C)(=C)", "C(-C)(-C)(=N)", "C(-C)(-C)(=O)", "C(-C)(-Cl)(=O)",
      "C(-C)(-H)(=C)", "C(-C)(-H)(=N)", "C(-C)(-H)(=O)", "C(-C)(-N)(=C)", "C(-C)(-N)(=N)", "C(-C)(-N)(=O)",
      "C(-C)(-O)(=O)", "C(-C)(=C)", "C(-C)(=N)", "C(-C)(=O)", "C(-Cl)(=O)", "C(-H)(-N)(=C)", "C(-H)(=C)", "C(-H)(=N)",
      "C(-H)(=O)", "C(-N)(=C)", "C(-N)(=N)", "C(-N)(=O)", "C(-O)(=O)", "N(-C)(=C)", "N(-C)(=O)", "N(-O)(=O)",
      "P(-O)(=O)", "S(-C)(=O)", "S(-O)(=O)", "S(=O)(=O)",
      ///Section 6/// // starts at bit 517
      "C-C-C#C", "O-C-C=N", "O-C-C=O", "N:C-S-H", "N-C-C=C", "O=S-C-C", "N#C-C=C", "C=N-N-C", "O=S-C-N", "S-S-C:C",
      "C:C-C=C", "S:C:C:C", "C:N:C-C", "S-C:N:C", "S:C:C:N", "S-C=N-C", "C-O-C=C", "N-N-C:C", "S-C=N-H", "S-C-S-C",
      "C:S:C-C", "O-S-C:C", "C:N-C:C", "N-S-C:C", "N-C:N:C", "N:C:C:N", "N-C:N:N", "N-C=N-C", "N-C=N-H", "N-C-S-C",
      "C-C-C=C", "C-N:C-H", "N-C:O:C", "O=C-C:C", "O=C-C:N", "C-N-C:C", "N:N-C-H", "O-C:C:N", "O-C=C-C", "N-C:C:N",
      "C-S-C:C", "Cl-C:C-C", "N-C=C-H", "Cl-C:C-H", "N:C:N-C", "Cl-C:C-O", "C-C:N:C", "C-C-S-C", "S=C-N-C", "Br-C:C-C",
      "H-N-N-H", "S=C-N-H", "C-O-H", "S:C:C-H",
      "O-N-C-C", "N-N-C-C", "H-C=C-H", "N-N-C-N", "O=C-N-N", "N=C-N-C", "C=C-C:C", "C:N-C-H", "C-N-N-H", "N:C:C-C",
      "C-C=C-C", "C:C-H", "Cl-C:C-Cl", "C:C:N-H", "H-N-C-H", "Cl-C-C-Cl", "N:C-C:C", "S-C:C-C", "S-C:C-H", "S-C:C-N",
      "S-C:C-O", "O=C-C-C", "O=C-C-N", "O=C-C-O", "N=C-C-C", "N=C-C-H", "C-N-C-H", "O-C:C-C", "O-C:C-H", "O-C:C-N",
      "O-C:C-O", "N-C:C-C", "N-C:C-H", "N-C:C-N", "O-C-C:C", "N-C-C:C", "Cl-C-C-C", "Cl-C-C-O", "C:C-C:C", "O=C-C=C",
      "Br-C-C-C", "N=C-C=C", "C=C-C-C", "N:C-O-H", "O=N-C:C", "O-C-N-H", "N-C-N-C", "Cl-C-C=O", "Br-C-C=O", "O-C-O-C",
      "C=C-C=C", "C:C-O-C", "O-C-C-N", "O-C-C-O", "N#C-C-C", "N-C-C-N", "C:C-C-C", "H-C-O-H", "N:C:N:C", "O-C-C=C",
      "O-C-C:C-C", "O-C-C:C-O", "N=C-C:C-H", "C:C-N-C:C", "C-C:C-C:C", "O=C-C-C-C", "O=C-C-C-N", "O=C-C-C-O",
      "C-C-C-C-C", "Cl-C:C-O-C", "C:C-C=C-C", "C-C:C-N-C", "C-S-C-C-C", "N-C:C-O-H", "O=C-C-C=O", "C-C:C-O-C",
      "C-C:C-O-H", "Cl-C-C-C-C", "N-C-C-C-C", "N-C-C-C-N", "C-O-C-C=C", "C:C-C-C-C", "N=C-N-C-C", "O=C-C-C:C",
      "Cl-C:C:C-C", "H-C-C=C-H", "N-C:C:C-C", "N-C:C:C-N", "O=C-C-N-C", "C-C:C:C-C", "C-O-C-C:C", "O=C-C-O-C",
      "O-C:C-C-C", "N-C-C-C:C", "C-C-C-C:C", "Cl-C-C-N-C", "C-O-C-O-C", "N-C-C-N-C", "N-C-O-C-C", "C-N-C-C-C",
      "C-C-O-C-C", "N-C-C-O-C", "C:C:N:N:C", "C-C-C-O-H", "C:C-C-C:C", "O-C-C=C-C", "C:C-O-C-C", "N-C:C:C:N",
      "O=C-O-C:C", "O=C-C:C-C", "O=C-C:C-N", "O=C-C:C-O", "C-O-C:C-C", "O=C:C:C", "C-N-C-C:C", "S-C:C:C-N", "O-C:C-O-C",
      "O-C:C-O-H", "C-C-O-C:C", "N-C-C:C-C", "C-C-C:C-C", "N-N-C-N-H", "C-N-C-N-C", "O-C-C-C-C", "O-C-C-C-N",
      "O-C-C-C-O", "C=C-C-C-C", "O-C-C-C=C", "O-C-C-C=O", "H-C-C-N-H", "C-C=N-N-C", "O=C-N-C-C", "O=C-N-C-H",
      "O=C-N-C-N", "O=N-C:C-N", "O=N-C:C-O", "O=C-N-C=O", "O-C:C:C-C", "O-C:C:C-N", "O-C:C:C-O", "N-C-N-C-C",
      "O-C-C-C:C", "C-C-N-C-C", "C-N-C:C-C", "C-C-S-C-C", "O-C-C-N-C", "C-C=C-C-C", "O-C-O-C-C", "O-C-C-O-C",
      "O-C-C-O-H", "C-C=C-C=C", "N-C:C-C-C", "C=C-C-O-C", "C=C-C-O-H", "C-C:C-C-C", "Cl-C:C-C=O", "Br-C:C:C-C",
      "O=C-C=C-C", "O=C-C=C-H", "O=C-C=C-N", "N-C-N-C:C", "Br-C-C-C:C", "N#C-C-C-C", "C-C=C-C:C", "C-C-C=C-C",
      "C-C-C-C-C-C", "O-C-C-C-C-C", "O-C-C-C-C-O", "O-C-C-C-C-N", "N-C-C-C-C-C", "O=C-C-C-C-C", "O=C-C-C-C-N",
      "O=C-C-C-C-O", "O=C-C-C-C=O", "C-C-C-C-C-C-C", "O-C-C-C-C-C-C", "O-C-C-C-C-C-O", "O-C-C-C-C-C-N", "O=C-C-C-C-C-C",
      "O=C-C-C-C-C-O", "O=C-C-C-C-C=O", "O=C-C-C-C-C-N", "C-C-C-C-C-C-C-C", "C-C-C-C-C-C(-C)-C", "O-C-C-C-C-C-C-C",
      "O-C-C-C-C-C(-C)-C", "O-C-C-C-C-C-O-C", "O-C-C-C-C-C(-O)-C", "O-C-C-C-C-C-N-C", "O-C-C-C-C-C(-N)-C",
      "O=C-C-C-C-C-C-C", "O=C-C-C-C-C(-O)-C", "O=C-C-C-C-C(=O)-C", "O=C-C-C-C-C(-N)-C", "C-C(-C)-C-C", "C-C(-C)-C-C-C",
      "C-C-C(-C)-C-C", "C-C(-C)(-C)-C-C", "C-C(-C)-C(-C)-C", "Cc1ccc(C)cc1", "Cc1ccc(O)cc1", "Cc1ccc(S)cc1",
      "Cc1ccc(N)cc1", "Cc1ccc(Cl)cc1", "Cc1ccc(Br)cc1", "Oc1ccc(O)cc1", "Oc1ccc(S)cc1", "Oc1ccc(N)cc1", "Oc1ccc(Cl)cc1",
      "Oc1ccc(Br)cc1", "Sc1ccc(S)cc1", "Sc1ccc(N)cc1", "Sc1ccc(Cl)cc1", "Sc1ccc(Br)cc1", "Nc1ccc(N)cc1",
      "Nc1ccc(Cl)cc1", "Nc1ccc(Br)cc1", "Clc1ccc(Cl)cc1", "Clc1ccc(Br)cc1", "Brc1ccc(Br)cc1", "Cc1cc(C)ccc1",
      "Cc1cc(O)ccc1", "Cc1cc(S)ccc1", "Cc1cc(N)ccc1", "Cc1cc(Cl)ccc1", "Cc1cc(Br)ccc1", "Oc1cc(O)ccc1", "Oc1cc(S)ccc1",
      "Oc1cc(N)ccc1", "Oc1cc(Cl)ccc1", "Oc1cc(Br)ccc1", "Sc1cc(S)ccc1", "Sc1cc(N)ccc1", "Sc1cc(Cl)ccc1",
      "Sc1cc(Br)ccc1", "Nc1cc(N)ccc1", "Nc1cc(Cl)ccc1", "Nc1cc(Br)ccc1", "Clc1cc(Cl)ccc1", "Clc1cc(Br)ccc1",
      "Brc1cc(Br)ccc1", "Cc1c(C)cccc1", "Cc1c(O)cccc1", "Cc1c(S)cccc1", "Cc1c(N)cccc1", "Cc1c(Cl)cccc1",
      "Cc1c(Br)cccc1", "Oc1c(O)cccc1", "Oc1c(S)cccc1", "Oc1c(N)cccc1", "Oc1c(Cl)cccc1", "Oc1c(Br)cccc1", "Sc1c(S)cccc1",
      "Sc1c(N)cccc1", "Sc1c(Cl)cccc1", "Sc1c(Br)cccc1", "Nc1c(N)cccc1", "Nc1c(Cl)cccc1", "Nc1c(Br)cccc1",
      "Clc1c(Cl)cccc1", "Clc1c(Br)cccc1", "Brc1c(Br)cccc1", "CC1CCC(C)CC1", "CC1CCC(O)CC1", "CC1CCC(S)CC1",
      "CC1CCC(N)CC1", "CC1CCC(Cl)CC1", "CC1CCC(Br)CC1", "OC1CCC(O)CC1", "OC1CCC(S)CC1", "OC1CCC(N)CC1", "OC1CCC(Cl)CC1",
      "OC1CCC(Br)CC1", "SC1CCC(S)CC1", "SC1CCC(N)CC1", "SC1CCC(Cl)CC1", "SC1CCC(Br)CC1", "NC1CCC(N)CC1",
      "NC1CCC(Cl)CC1", "NC1CCC(Br)CC1", "ClC1CCC(Cl)CC1", "ClC1CCC(Br)CC1", "BrC1CCC(Br)CC1", "CC1CC(C)CCC1",
      "CC1CC(O)CCC1", "CC1CC(S)CCC1", "CC1CC(N)CCC1", "CC1CC(Cl)CCC1", "CC1CC(Br)CCC1", "OC1CC(O)CCC1", "OC1CC(S)CCC1",
      "OC1CC(N)CCC1", "OC1CC(Cl)CCC1", "OC1CC(Br)CCC1", "SC1CC(S)CCC1", "SC1CC(N)CCC1", "SC1CC(Cl)CCC1",
      "SC1CC(Br)CCC1", "NC1CC(N)CCC1", "NC1CC(Cl)CCC1", "NC1CC(Br)CCC1", "ClC1CC(Cl)CCC1", "ClC1CC(Br)CCC1",
      "BrC1CC(Br)CCC1", "CC1C(C)CCCC1", "CC1C(O)CCCC1", "CC1C(S)CCCC1", "CC1C(N)CCCC1", "CC1C(Cl)CCCC1",
      "CC1C(Br)CCCC1", "OC1C(O)CCCC1", "OC1C(S)CCCC1", "OC1C(N)CCCC1", "OC1C(Cl)CCCC1", "OC1C(Br)CCCC1", "SC1C(S)CCCC1",
      "SC1C(N)CCCC1", "SC1C(Cl)CCCC1", "SC1C(Br)CCCC1", "NC1C(N)CCCC1", "NC1C(Cl)CCCC1", "NC1C(Br)CCCC1",
      "ClC1C(Cl)CCCC1", "ClC1C(Br)CCCC1", "BrC1C(Br)CCCC1", "CC1CC(C)CC1", "CC1CC(O)CC1", "CC1CC(S)CC1", "CC1CC(N)CC1",
      "CC1CC(Cl)CC1", "CC1CC(Br)CC1", "OC1CC(O)CC1", "OC1CC(S)CC1", "OC1CC(N)CC1", "OC1CC(Cl)CC1", "OC1CC(Br)CC1",
      "SC1CC(S)CC1", "SC1CC(N)CC1", "SC1CC(Cl)CC1", "SC1CC(Br)CC1", "NC1CC(N)CC1", "NC1CC(Cl)CC1", "NC1CC(Br)CC1",
      "ClC1CC(Cl)CC1", "ClC1CC(Br)CC1", "BrC1CC(Br)CC1", "CC1C(C)CCC1", "CC1C(O)CCC1", "CC1C(S)CCC1", "CC1C(N)CCC1",
      "CC1C(Cl)CCC1", "CC1C(Br)CCC1", "OC1C(O)CCC1", "OC1C(S)CCC1", "OC1C(N)CCC1", "OC1C(Cl)CCC1", "OC1C(Br)CCC1",
      "SC1C(S)CCC1", "SC1C(N)CCC1", "SC1C(Cl)CCC1", "SC1C(Br)CCC1", "NC1C(N)CCC1", "NC1C(Cl)CC1", "NC1C(Br)CCC1",
      "ClC1C(Cl)CCC1", "ClC1C(Br)CCC1", "BrC1C(Br)CCC1"
    };

    // atom number atom symbol; numbers come from PubChem Substructure Fingerprint Description
    // numbers represent the log2 max value of interest for the element type
    const size_t MACCS::s_ElementArray[ 92] =
    {
      // H                                                  He
      4,                                                 0,
      // Li Be                               B  C  N  O  F  Ne
      2, 1,                               3, 5, 4, 5, 3, 0,
      // Na Mg                               Al Si P  S  Cl Ar
      2, 1,                               1, 2, 3, 4, 4, 0,
      // K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
      2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 0,
      // Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 0,
      // Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
      // Fr Ra Ac Th Pa U
      0, 0, 0, 0, 0, 1
    };
  } // namespace descriptor
} // namespace bcl
