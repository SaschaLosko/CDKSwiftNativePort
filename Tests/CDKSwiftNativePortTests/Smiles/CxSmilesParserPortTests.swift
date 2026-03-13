import XCTest
@testable import CDKSwiftNativePort

final class CxSmilesParserPortTests: XCTestCase {

    func testParsesCxAtomLabelsLayer() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C* |$;R1$|")

        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.atoms[0].element, "C")
        XCTAssertEqual(molecule.atoms[1].element, "R1")
    }

    func testParsesCxRacemicLayerWithoutFailure() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC[C@H](O)[C@H](O)CCCCCC |r|")

        XCTAssertEqual(molecule.atomCount, 12)
        XCTAssertTrue(molecule.atoms.contains(where: { $0.chirality != .none }))
    }

    func testParsesCxRacemicFragmentIndicesLayer() throws {
        let split = try CDKCxSmilesParser.split("CC[C@H](O)[C@H](O)CCCCCC |r:1|", enabled: true)
        XCTAssertEqual(split.state.racemicFragments, [1])
    }

    func testParsesCxFragmentLayerWithoutFailure() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C.*.C |f:0.2|")
        XCTAssertEqual(molecule.atomCount, 3)
    }

    func testParsesCxWithTitle() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |$C;C$| ethane")
        XCTAssertEqual(molecule.name, "ethane")
    }

    func testRejectsMalformedCxLayer() {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        XCTAssertThrowsError(try parser.parseSmiles("CC |$;R1$"))
    }

    func testRejectsMalformedCxRacemicFragmentLayer() {
        XCTAssertThrowsError(try CDKCxSmilesParser.split("CC |r:a|", enabled: true))
    }

    func testParsesMarkushRGroupsFromUserExample() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")

        XCTAssertEqual(connectedComponents(in: molecule).count, 4)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == "R1" && $0.attachmentPoint == 1 }.count, 3)
        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups.first?.subscriptText, "1-3")
    }

    func testParsesNestedCxInRGroupDefinitions() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C* |$;R1$,RG:_R1={Cl* |$;_AP1$|},{Br* |$;_AP1$|}|")

        XCTAssertEqual(connectedComponents(in: molecule).count, 3)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == "R1" && $0.attachmentPoint == 1 }.count, 2)
        XCTAssertTrue(molecule.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Cl" })
        XCTAssertTrue(molecule.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Br" })
    }

    func testParsesSimpleLinkNodeIntoRepeatUnitAnnotation() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3|")

        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups[0].kind, .structureRepeatUnit)
        XCTAssertEqual(molecule.sgroups[0].subscriptText, "1-3")
        XCTAssertEqual(molecule.sgroups[0].atomIDs.count, 1)
    }

    func testParsesMultipleLinkNodes() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3,3:1.4|")

        XCTAssertEqual(molecule.sgroups.count, 2)
        XCTAssertEqual(Set(molecule.sgroups.map(\.subscriptText)), Set(["1-3", "1-4"]))
    }

    func testParsesLinkNodeWithAttachedAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

        let molecule = try parser.parseSmiles("C1NC(O)NC1 |LN:2:1.3|")
        XCTAssertEqual(molecule.sgroups.count, 1)
        let linkedElements = molecule.sgroups[0].atomIDs.compactMap { atomID in
            molecule.atom(id: atomID)?.element
        }
        XCTAssertEqual(Set(linkedElements), Set(["C", "O"]))

        let withExplicitCrossingBonds = try parser.parseSmiles("C1NC(O)NC1 |LN:2:1.3.1.3|")
        XCTAssertEqual(withExplicitCrossingBonds.sgroups.count, 1)
        XCTAssertEqual(withExplicitCrossingBonds.sgroups[0].atomIDs.count, 2)
        XCTAssertEqual(Set(withExplicitCrossingBonds.sgroups[0].crossingBondIDs).count, 2)
    }

    func testPreservesDisconnectedMarkushAlternativesUsingComponentGroupIDs() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("*c1ncccc1 |$R1$,RG:_R1={Cl},{Br},{C(=O)*.Cl |$;;R2$,RG:_R2={OC},{N}|}|")

        let r1Atoms = molecule.atoms.filter { $0.rGroupMembership == "R1" }
        let grouped = Dictionary(grouping: r1Atoms, by: { $0.componentGroupID ?? -1 })

        XCTAssertGreaterThanOrEqual(grouped.count, 3)
        XCTAssertTrue(grouped.values.contains { atoms in
            atoms.contains(where: { $0.element == "Cl" }) &&
            atoms.contains(where: { $0.element == "C" || $0.element == "O" })
        })
    }

    private func connectedComponents(in molecule: Molecule) -> [Set<Int>] {
        var adjacency: [Int: Set<Int>] = [:]
        for atom in molecule.atoms {
            _ = adjacency[atom.id, default: []]
        }
        for bond in molecule.bonds {
            adjacency[bond.a1, default: []].insert(bond.a2)
            adjacency[bond.a2, default: []].insert(bond.a1)
        }

        var seen = Set<Int>()
        var components: [Set<Int>] = []

        for atom in molecule.atoms where !seen.contains(atom.id) {
            var queue = [atom.id]
            var component = Set<Int>()
            seen.insert(atom.id)

            while let current = queue.popLast() {
                component.insert(current)
                for neighbor in adjacency[current, default: []] where !seen.contains(neighbor) {
                    seen.insert(neighbor)
                    queue.append(neighbor)
                }
            }

            components.append(component)
        }

        return components
    }
}
