import XCTest
@testable import CDKSwiftNativePort

final class MDLV3000ReaderScaffoldTests: XCTestCase {

    func testEmptyStringThrowsExpectedHeaderError() {
        XCTAssertThrowsError(try CDKMDLV3000Reader.read(text: "")) { error in
            guard case let ChemError.parseFailed(message) = error else {
                return XCTFail("Unexpected error: \(error)")
            }
            XCTAssertEqual(message, "Expected a header line, but found nothing.")
        }
    }

    func testParsesV3000Scaffold() throws {
        let scaffold = try CDKMDLV3000Reader.parseScaffold(text: minimalV3000)

        XCTAssertEqual(scaffold.title, "V3000Example")
        XCTAssertEqual(scaffold.atomCount, 2)
        XCTAssertEqual(scaffold.bondCount, 1)
        XCTAssertEqual(scaffold.atomRecords.count, 2)
        XCTAssertEqual(scaffold.bondRecords.count, 1)
        XCTAssertEqual(scaffold.atomRecords[0].symbol, "C")
        XCTAssertEqual(scaffold.bondRecords[0].type, "1")
        XCTAssertEqual(scaffold.bondRecords[0].a1, 1)
        XCTAssertEqual(scaffold.bondRecords[0].a2, 2)
    }

    func testReadParsesBasicV3000IntoMolecule() throws {
        let molecule = try CDKMDLV3000Reader.read(text: minimalV3000)

        XCTAssertEqual(molecule.name, "V3000Example")
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.atoms.map(\.element), ["C", "O"])
        XCTAssertEqual(molecule.bonds[0].order, .single)
        XCTAssertEqual(molecule.bonds[0].a1, 1)
        XCTAssertEqual(molecule.bonds[0].a2, 2)
    }

    func testMDLReaderDispatchesV3000ToMolecule() throws {
        let molecule = try CDKMDLReader.read(text: minimalV3000)
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
    }

    func testParsesV3000AtomAndBondAttributes() throws {
        let molecule = try CDKMDLV3000Reader.read(text: attributedV3000)
        XCTAssertEqual(molecule.atomCount, 4)
        XCTAssertEqual(molecule.bondCount, 3)

        let atom1 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 1 }))
        XCTAssertEqual(atom1.charge, 1)
        XCTAssertEqual(atom1.isotopeMassNumber, 13)
        XCTAssertEqual(atom1.explicitHydrogenCount, 0)

        let atom2 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 2 }))
        XCTAssertEqual(atom2.radical, 2)

        let atom3 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 3 }))
        XCTAssertEqual(atom3.element, "L")
        XCTAssertEqual(atom3.queryType, .anyAtom)
        XCTAssertEqual(atom3.atomList ?? [], ["N", "O"])

        let atom4 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 4 }))
        XCTAssertEqual(atom4.queryType, .anyNonHydrogen)

        let bond1 = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 1 }))
        XCTAssertEqual(bond1.queryType, .singleOrAromatic)
        XCTAssertEqual(bond1.stereo, .up)

        let bond2 = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 2 }))
        XCTAssertEqual(bond2.queryType, .doubleOrAromatic)

        let bond3 = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 3 }))
        XCTAssertEqual(bond3.queryType, .any)
    }

    func testParsesContinuationLines() throws {
        let molecule = try CDKMDLV3000Reader.read(text: continuationV3000)
        let atom = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 1 }))

        XCTAssertEqual(atom.element, "L")
        XCTAssertEqual(atom.atomList ?? [], ["C", "N", "O"])
        XCTAssertEqual(atom.queryType, .anyAtom)
    }

    func testGeneratesLayoutWhenCoordinatesAreMissing() throws {
        let molecule = try CDKMDLV3000Reader.read(text: zeroCoordinateV3000)
        let box = try XCTUnwrap(molecule.boundingBox())
        XCTAssertGreaterThan(max(box.width, box.height), 0.1)
    }

    func testParseScaffoldRejectsMissingCTAB() {
        XCTAssertThrowsError(try CDKMDLV3000Reader.parseScaffold(text: invalidV3000))
    }

    func testAppliesCollectionAndSGroupSemantics() throws {
        let scaffold = try CDKMDLV3000Reader.parseScaffold(text: collectionAndSGroupV3000)
        XCTAssertEqual(scaffold.collectionLines.count, 1)
        XCTAssertEqual(scaffold.sgroupLines.count, 2)

        let molecule = try CDKMDLV3000Reader.read(text: collectionAndSGroupV3000)
        let bond = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 1 }))
        XCTAssertEqual(bond.order, .double)
        XCTAssertEqual(bond.stereo, .either)

        let atom = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 1 }))
        XCTAssertEqual(atom.element, "Me")
        XCTAssertEqual(atom.attachmentPoint, 1)
        XCTAssertEqual(molecule.name, "V3000CollectionAndSGroup")
        XCTAssertEqual(molecule.orderedDataFieldNames, ["CDK"])
        XCTAssertEqual(molecule.dataFieldValues(named: "CDK"), ["AliasName"])
    }

    func testParsesTrailingSDFFieldsFromV3000Record() throws {
        let molecule = try CDKMDLV3000Reader.read(text: annotatedV3000)

        XCTAssertEqual(molecule.orderedDataFieldNames, ["ID", "Tags"])
        XCTAssertEqual(molecule.dataFieldValues(named: "ID"), ["001"])
        XCTAssertEqual(molecule.dataFieldValues(named: "Tags"), ["alpha", "beta"])
    }

    func testParsesTrailingSDFFieldsWithoutDelimiterAndWithDelimiterNoise() throws {
        let noDelimiter = try CDKMDLV3000Reader.read(text: annotatedV3000NoDelimiter)
        XCTAssertEqual(noDelimiter.orderedDataFieldNames, ["ID", "Tags"])
        XCTAssertEqual(noDelimiter.dataFieldValues(named: "Tags"), ["alpha", "beta"])

        let noisyDelimiter = try CDKMDLV3000Reader.read(text: annotatedV3000DelimiterNoise)
        XCTAssertEqual(noisyDelimiter.orderedDataFieldNames, ["ID", "Tags"])
        XCTAssertEqual(noisyDelimiter.dataFieldValues(named: "Tags"), ["alpha", "beta"])
    }

    func testParsesHILITECollectionIntoMoleculeSelectionState() throws {
        let molecule = try CDKMDLV3000Reader.read(text: highlightedV3000)

        XCTAssertEqual(molecule.highlightedAtomIDs, [2, 3, 5, 4])
        XCTAssertEqual(molecule.highlightedBondIDs, [2])
    }

    func testParsesBondTopologyAndQueryTypes() throws {
        let molecule = try CDKMDLV3000Reader.read(text: bondTopologyV3000)

        XCTAssertEqual(molecule.bondCount, 4)
        XCTAssertEqual(molecule.bonds[0].queryType, .singleOrDouble)
        XCTAssertEqual(molecule.bonds[0].topology, .ring)
        XCTAssertEqual(molecule.bonds[1].queryType, .doubleOrAromatic)
        XCTAssertEqual(molecule.bonds[1].topology, .chain)
        XCTAssertEqual(molecule.bonds[2].queryType, .any)
        XCTAssertNil(molecule.bonds[2].topology)
        XCTAssertEqual(molecule.bonds[3].order, .single)
        XCTAssertEqual(molecule.bonds[3].topology, .ring)
    }

    func testIgnoresInvalidBondTopologyValue() throws {
        let molecule = try CDKMDLV3000Reader.read(text: invalidBondTopologyV3000)

        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.bonds[0].queryType, .singleOrDouble)
        XCTAssertNil(molecule.bonds[0].topology)
    }

    func testRejectsUnsupportedAndInvalidBondTypes() {
        XCTAssertThrowsError(try CDKMDLV3000Reader.read(text: unsupportedBondType9V3000))
        XCTAssertThrowsError(try CDKMDLV3000Reader.read(text: unsupportedBondType10V3000))
        XCTAssertThrowsError(try CDKMDLV3000Reader.read(text: invalidBondType11V3000))
    }

    func testParsesNonSequentialAtomIndices() throws {
        let molecule = try CDKMDLV3000Reader.read(text: nonSequentialAtomIndexV3000)

        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.atoms.map(\.id), [1, 22])
        XCTAssertEqual(molecule.bonds.first?.a1, 1)
        XCTAssertEqual(molecule.bonds.first?.a2, 22)
    }

    func testParsesAtomMappingAndHydrogenIsotopes() throws {
        let molecule = try CDKMDLV3000Reader.read(text: mappedIsotopeV3000)

        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(molecule.atoms[0].atomMapNumber, 7)
        XCTAssertEqual(molecule.atoms[0].atomClass, 7)
        XCTAssertEqual(molecule.atoms[1].element, "H")
        XCTAssertEqual(molecule.atoms[1].isotopeMassNumber, 2)
        XCTAssertEqual(molecule.atoms[2].element, "H")
        XCTAssertEqual(molecule.atoms[2].isotopeMassNumber, 3)
    }

    func testParsesPseudoAtomLabelsAndPreservesBondReferences() throws {
        let molecule = try CDKMDLV3000Reader.read(text: pseudoAtomV3000)

        XCTAssertEqual(molecule.atomCount, 10)
        let pseudo = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 10 }))
        XCTAssertEqual(pseudo.element, "Leu")
        XCTAssertTrue(molecule.bonds.contains(where: { $0.a1 == 7 && $0.a2 == 10 }))
    }

    func testParsesAromaticBondTypeFourAsAromaticBondOrder() throws {
        let molecule = try CDKMDLV3000Reader.read(text: aromaticBondTypeFourV3000)

        XCTAssertEqual(molecule.atomCount, 6)
        XCTAssertEqual(molecule.bondCount, 6)
        XCTAssertTrue(molecule.bonds.allSatisfy { $0.order == .aromatic })
    }

    func testParsesZeroDimStereoAtomsFromCFGAttributes() throws {
        let molecule = try CDKMDLV3000Reader.read(text: stereo0dV3000)
        let chiralAtoms = molecule.atoms.filter { $0.chirality != .none }

        XCTAssertEqual(chiralAtoms.map(\.id), [5, 6])
        XCTAssertTrue(chiralAtoms.allSatisfy {
            $0.cxStereoGroup == CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 1)
        })
    }

    func testParsesStereoCollectionsAndChiralFlag() throws {
        let absMolecule = try CDKMDLV3000Reader.read(text: chiralFlagV3000)
        XCTAssertEqual(absMolecule.atoms.first?.cxStereoGroup,
                       CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0))

        let mixedMolecule = try CDKMDLV3000Reader.read(text: stereoCollectionsV3000)
        XCTAssertEqual(mixedMolecule.atoms[0].cxStereoGroup,
                       CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0))
        XCTAssertEqual(mixedMolecule.atoms[1].cxStereoGroup,
                       CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 1))
        XCTAssertEqual(mixedMolecule.atoms[2].cxStereoGroup,
                       CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 1))
    }

    func testParsesPositionVariationBondAttributesIntoExtMulticenterSgroup() throws {
        let molecule = try CDKMDLV3000Reader.read(text: positionalVariationV3000)

        let multicenter = try XCTUnwrap(molecule.sgroups.first(where: { $0.kind == .extMulticenter }))
        XCTAssertEqual(multicenter.crossingBondIDs, [8])
        XCTAssertEqual(multicenter.atomIDs, [8, 2, 3, 4, 5, 6])
    }

    func testParsesRepeatedBracketAndParentSgroupAttributes() throws {
        let molecule = try CDKMDLV3000Reader.read(text: nestedSgroupsV3000)

        XCTAssertEqual(molecule.sgroups.count, 2)
        let parent = molecule.sgroups[0]
        let child = molecule.sgroups[1]
        XCTAssertEqual(parent.keyword, "FOR")
        XCTAssertEqual(parent.brackets.count, 2)
        XCTAssertEqual(parent.childGroupIndices, [1])
        XCTAssertEqual(child.componentNumber, 1)
        XCTAssertEqual(child.parentAtomIDs, [2])
    }

    private let minimalV3000 = """
V3000Example
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 O 1.2000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let attributedV3000 = """
V3000Attributes
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0 CHG=1 MASS=13 HCOUNT=0
M  V30 2 O 1.2000 0.0000 0.0000 0 RAD=2
M  V30 3 [N,O] 2.4000 0.0000 0.0000 0
M  V30 4 A 3.6000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 6 1 2 CFG=1
M  V30 2 7 2 3
M  V30 3 8 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let continuationV3000 = """
V3000Continuation
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 L 0.0000 0.0000 0.0000 0 ATOMLIST=(3 C N -
M  V30 O)
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let zeroCoordinateV3000 = """
V3000ZeroCoords
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let invalidV3000 = """
V3000Invalid
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  END
"""

    private let collectionAndSGroupV3000 = """
V3000CollectionAndSGroup
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 C 1.3000 0.0000 0.0000 0
M  V30 3 O 2.6000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS BONDS=(1 1)
M  V30 END COLLECTION
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(1 1) LABEL="Me" ATTCHPT=1
M  V30 2 DAT 0 ATOMS=(1 3) FIELDNAME="CDK" FIELDDATA="AliasName"
M  V30 END SGROUP
M  V30 END CTAB
M  END
"""

    private let annotatedV3000 = """
V3000Annotated
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
>  <ID>
001

>  <Tags>
alpha
beta
"""

    private let annotatedV3000NoDelimiter = """
V3000Annotated
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
>  <ID>
001

>  <Tags>
alpha
beta
"""

    private let annotatedV3000DelimiterNoise = """
V3000Annotated
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
>  <ID>
001

>  <Tags>
alpha
beta

$$$$    <--- if we see that, it's enough and anything here is ignored
"""

    private let highlightedV3000 = """
V3000Highlighted
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.5000 0.0000 0.0000 0
M  V30 2 C 0.7500 1.2990 0.0000 0
M  V30 3 C -0.7500 1.2990 0.0000 0
M  V30 4 C -1.5000 0.0000 0.0000 0
M  V30 5 C -0.7500 -1.2990 0.0000 0
M  V30 6 C 0.7500 -1.2990 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 6 1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/HILITE ATOMS=(4 2 3 5 4) BONDS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
"""

    private let bondTopologyV3000 = """
V3000BondTopology
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 3 C 2 0 0 0
M  V30 4 C 3 0 0 0
M  V30 5 C 4 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 5 1 2 TOPO=1
M  V30 2 7 2 3 TOPO=2
M  V30 3 8 3 4
M  V30 4 1 4 5 TOPO=1
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let invalidBondTopologyV3000 = """
V3000InvalidTopology
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 N 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 5 1 2 TOPO=3
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let unsupportedBondType9V3000 = """
V3000Unsupported9
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 9 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let unsupportedBondType10V3000 = """
V3000Unsupported10
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 10 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let invalidBondType11V3000 = """
V3000Invalid11
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 11 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let nonSequentialAtomIndexV3000 = """
V3000NonSequential
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 22 C 1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 22
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let mappedIsotopeV3000 = """
V3000MappedIsotope
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 7
M  V30 2 D 1 0 0 0
M  V30 3 T 2 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let pseudoAtomV3000 = """

  Marvin  01211213222D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.925 0.935 0 0
M  V30 2 C -3.2587 0.165 0 0
M  V30 3 C -3.2587 -1.375 0 0
M  V30 4 C -1.925 -2.145 0 0
M  V30 5 C -0.5913 -1.375 0 0
M  V30 6 C -0.5913 0.165 0 0
M  V30 7 C 0.8962 0.5636 0 0
M  V30 8 C 2.2299 -0.2064 0 0
M  V30 9 C 3.7174 0.1922 0 0
M  V30 10 Leu 0.8962 2.1036 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 1 6
M  V30 3 2 2 3
M  V30 4 1 3 4
M  V30 5 2 4 5
M  V30 6 1 5 6
M  V30 7 1 6 7
M  V30 8 1 7 8
M  V30 9 1 8 9
M  V30 10 1 7 10
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let aromaticBondTypeFourV3000 = """

  Mrv2403 08072416472D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.5415 1.04 0 0
M  V30 2 C -5.8749 0.2701 0 0
M  V30 3 C -5.8749 -1.27 0 0
M  V30 4 C -4.5415 -2.04 0 0
M  V30 5 C -3.2078 -1.27 0 0
M  V30 6 C -3.2078 0.2701 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 4 5
M  V30 5 4 5 6
M  V30 6 4 6 1
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let chiralFlagV3000 = """
V3000ChiralFlag
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0 CFG=2
M  V30 2 O 1 0 0 0
M  V30 3 N -1 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let stereo0dV3000 = """

  ChemDraw0207221100

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 N 0.0 0.0 0.0 0
M  V30 3 C 0.0 0.0 0.0 0
M  V30 4 C 0.0 0.0 0.0 0
M  V30 5 C 0.0 0.0 0.0 0 CFG=1
M  V30 6 C 0.0 0.0 0.0 0 CFG=1
M  V30 7 O 0.0 0.0 0.0 0
M  V30 8 O 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 1 5 8
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let stereoCollectionsV3000 = """
V3000StereoCollections
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0 CFG=2
M  V30 2 C 1 0 0 0 CFG=1
M  V30 3 C 2 0 0 0 CFG=2
M  V30 4 O 3 0 0 0
M  V30 5 O 4 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 4 CFG=1
M  V30 2 1 2 4 CFG=1
M  V30 3 1 3 5 CFG=1
M  V30 4 1 1 2
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 1)
M  V30 MDLV30/STERAC1 ATOMS=(1 2)
M  V30 MDLV30/STEREL1 ATOMS=(1 3)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
"""

    private let positionalVariationV3000 = """
V3000PositionalVariation
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1 0 0 0
M  V30 2 C 2 0 0 0
M  V30 3 C 3 0 0 0
M  V30 4 C 4 0 0 0
M  V30 5 C 5 0 0 0
M  V30 6 C 6 0 0 0
M  V30 7 C 7 0 0 0
M  V30 8 C 8 0 0 0
M  V30 9 * 9 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9 ATTACH=ANY ENDPTS=(5 2 3 4 5 6)
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let nestedSgroupsV3000 = """
V3000NestedSgroups
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 3 O 2 0 0 0
M  V30 4 C 3 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 FOR 0 ATOMS=(4 1 2 3 4) BRKXYZ=(9 -1 -1 0 -1 1 0 0 0 0) BRKXYZ=(9 4 1 0 4 -1 0 0 0 0)
M  V30 2 COM 0 ATOMS=(2 2 3) PARENT=1 PATOMS=(1 2) COMPNO=1
M  V30 END SGROUP
M  V30 END CTAB
M  END
"""
}
