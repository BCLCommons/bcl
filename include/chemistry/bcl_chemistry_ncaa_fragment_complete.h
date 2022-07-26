// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users: 
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users: 
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University: 
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the 
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University. 
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org 
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

#ifndef BCL_CHEMISTRY_NCAA_FRAGMENT_COMPLETE_H_
#define BCL_CHEMISTRY_NCAA_FRAGMENT_COMPLETE_H_

// include the namespace header
#include <biol/bcl_biol_atom_type_data.h>
#include <chemistry/bcl_chemistry_fragment_add_med_chem.h>
#include <function/bcl_function_member_unary_const.h>
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NCAAFragmentComplete
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_ncaa_fragment_complete.cpp @endlink
    //! @author vuot2, brownpb1, mendenjl
    //! @date Mar 14, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NCAAFragmentComplete :
      public FragmentComplete
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private :

      // automatically find the type of the backbone (only apply
      // if the backbone of the input NCAA structure is complete
      bool m_AutoFindBackboneType;
      // storage::Map< size_t, biol::AtomType> m_Backbone;    //!< Biol atom type and indices of the
      //storage::Map< size_t, biol::AtomType> m_FirstThreeSidechainAtoms;
      //! to append functional groups to backbones
      chemistry::FragmentAddMedChem m_AddMedChem;

      // properties of the NCAA
      storage::Map< std::string, std::string> m_Properties;

    public :
      // Different choices of backbone types
      enum Backbone_Type
      {
        e_Alpha,
        e_Peptoid,
        e_Nmethyl_Alpha,
        s_NumberBackboneTypes
      };

      //! @brief atom type as string
      //! @param ATOM_TYPE the name of the atom type
      //! @return the string for the atom type
      static const std::string &GetBackboneTypeName( const Backbone_Type &BACKBONE_TYPE);

      //! @brief Initialization enum I/O helper
      typedef util::WrapperEnum< Backbone_Type, &GetBackboneTypeName, s_NumberBackboneTypes> BackboneTypeEnum;

    private:

      typedef storage::Pair< const FragmentComplete&, const std::string> t_input;             //!> the wrapper parameter
      typedef storage::Pair< bool, chemistry::FragmentComplete > t_output;
      typedef storage::Pair
      <
          // Function that returns the dipeptide from input NCAA structure
          function::MemberUnaryConst< NCAAFragmentComplete, const t_input&, t_output>,
          // Function that returns the backbone properties for the dipeptide form
          function::MemberUnaryConst< NCAAFragmentComplete, const size_t &, storage::Map< std::string, std::string>>
      >
      t_functions;

      BackboneTypeEnum m_BackboneType;                                          //!> The type of atom type hash

    public:
      //! @brief Get the name, description and functions of the given backbone type
      //! @param BACKBONE_TYPE the atom type used to build the dipeptide form of the NCAA
      //! @return the short name or abbreviation of the class
      static const storage::Triplet< std::string, std::string, t_functions> &GetBackboneTypeInfo(
        const Backbone_Type &BACKBONE_TYPE);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      NCAAFragmentComplete();

      //! @brief Construct a NCAA from input a input NCAA structure
      //!        Automatically determines the type of the backbone based on the input NCAA
      //!        Should only be used when the backbone of the input NCAA are complete
      NCAAFragmentComplete
      (
        const FragmentComplete& INPUT_NCAA
      );

      //! @brief Construct a NCAA from input a input NCAA structure and a give type of the backbone
      NCAAFragmentComplete
      (
        const FragmentComplete& INPUT_NCAA,
        const Backbone_Type BACKBONE_TYPE
      );

      //! @brief pass-through constructor
      //! @param CONSTRUCTOR_OBJECT object that can be used to construct the base class
      template< typename t_DataType>
      NCAAFragmentComplete( const t_DataType &CONSTRUCTOR_OBJECT) :
        FragmentComplete( CONSTRUCTOR_OBJECT)
      {
        // note that we must call the conformation interfaces' implementation of GetNumberValences, since this class
        // overrides GetNumberValences, always returning 0; but here we need to check that this is true!
        // This code is taken from the molecule complete class, the NCAA might
        // not have valence = 0 if the NCAA is in the residue form
        /*if( ConformationInterface::GetNumberValences() != size_t( 0))
        {
          SaturateWithH();
          BCL_Assert( ConformationInterface::GetNumberValences() == size_t( 0), " Molecule should have 0 valences");
        }*/
      }

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      NCAAFragmentComplete *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the number of valences
      //! @return the number of valences (=0)
      size_t GetNumberValences() const
      {
        return 0;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // helpers   //
    ///////////////

    private :

      //! @brief automatically finds the type of backbone for the input NCAA structure
      //!        Should only be used in the case the backbone of the input NCAA is complete
      //! @return the chi1 atom index

      const Backbone_Type FindBackBoneType() const;

      //! @brief return the chi1 atom indices
      //! @param NCAA: the atom vector of NCAA
      //! @param BB_CONN_INDICES: the indices of three bb atoms that are closest to Chi 1 Atom
      //! @return the chi1 atom index
      const size_t FindChi1Index
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
        storage::Vector< size_t> BB_CONN_INDICES
      ) const;

      ///////////////
      // Alpha AAs //
      ///////////////

      //! @brief load neutral glycine residue from library
      //! @return the neutral glycine as the ncaa base
      const storage::Pair< bool, chemistry::FragmentComplete> ReadGlycineBase() const;

      //! @brief load neutral glycine dipeptide from library as backbone for alpha AA
      //! @return the neutral glycine dipeptide as the ncaa base
      const storage::Triplet< bool, size_t, chemistry::FragmentComplete> ReadGlycineDipeptideBackbone() const;

      //! @brief Find the CA chirality for alpha + Nmethyl alpha NCAAs
      //! @param NCAA: the atom vector of NCAA
      //! @param C_INDEX: the index of backbone C atom
      //! @param N_INDEX: the index of the backbone N atom
      //! @param CHI1_INDEX: the index of the chi1 angle atom
      const std::string FindAlphaCAChirarity
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
        const size_t &CA_INDEX,
        const size_t &C_INDEX,
        const size_t &N_INDEX,
        const size_t &CHI1_INDEX
      ) const;

      ///////////////
      // peptoid   //
      ///////////////

      //! @brief Find the CA chirality for alpha + Nmethyl alpha NCAAs
      //! @param NCAA: the atom vector of NCAA
      //! @param C_INDEX: the index of backbone C atom
      //! @param N_INDEX: the index of the backbone N atom
      //! @param CHI1_INDEX: the index of the chi1 angle atom
      const std::string FindPeptoidNChirarity
      (
        const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
        const size_t &C_INDEX,
        const size_t &N_INDEX,
        const size_t &CHI1_INDEX
      ) const;

      ////////////////////
      // N methyl alpha //
      ////////////////////

      //! @brief load neutral glycine residue from library
      //! @return the neutral glycine as the ncaa base
      const storage::Pair< bool, chemistry::FragmentComplete> ReadNMethylGlycineBase() const;

      //! @brief load neutral glycine dipeptide from library as backbone for alpha AA
      //! @return the neutral glycine dipeptide as the ncaa base
      const storage::Triplet< bool, size_t, chemistry::FragmentComplete> ReadNMethylGlycineDipeptideBackbone() const;


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

    }; // class NCAAFragmentComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_NCAA_FRAGMENT_COMPLETE_H_ 
