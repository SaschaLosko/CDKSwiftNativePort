import Foundation
import XCTest

final class PackageBoundaryTests: XCTestCase {
    func testSourcesContainNoAppLevelCouplingMarkers() throws {
        let root = try packageRoot()
        let sourcesRoot = root.appendingPathComponent("Sources/CDKSwiftNativePort", isDirectory: true)

        let forbiddenMarkers = [
            "de.losko.atomlens",
            "com_chemsketcher",
            "AtomLens",
            "import AtomLens",
            "import AtomLensQuickLookSupport",
            "import AtomLensSpotlightImporter",
            "import AppKit",
            "import UIKit",
            "import CoreSpotlight",
            "import QuickLook",
            "import QuickLookThumbnailing"
        ]

        var violations: [String] = []
        let enumerator = FileManager.default.enumerator(at: sourcesRoot,
                                                        includingPropertiesForKeys: nil)

        while let next = enumerator?.nextObject() as? URL {
            guard next.pathExtension == "swift" else { continue }
            let content = try String(contentsOf: next, encoding: .utf8)
            for marker in forbiddenMarkers where content.contains(marker) {
                let relativePath = next.path.replacingOccurrences(of: sourcesRoot.path + "/", with: "")
                violations.append("\(relativePath): contains '\(marker)'")
            }
        }

        XCTAssertTrue(violations.isEmpty, """
        Found package boundary violations:
        \(violations.joined(separator: "\n"))
        """)
    }

    func testPackageManifestDoesNotDependOnAtomLensTargets() throws {
        let root = try packageRoot()
        let manifestURL = root.appendingPathComponent("Package.swift")
        let manifest = try String(contentsOf: manifestURL, encoding: .utf8)

        let forbiddenMarkers = [
            "AtomLens",
            "AtomLensQuickLookSupport",
            "AtomLensSpotlightImporter",
            "path: \"../AtomLens",
            "path:\"../AtomLens"
        ]

        let violations = forbiddenMarkers.filter { manifest.contains($0) }
        XCTAssertTrue(violations.isEmpty, """
        Package manifest contains forbidden app-coupling markers:
        \(violations.joined(separator: "\n"))
        """)
    }

    func testWorkspaceChemistryLayerContainsOnlyAliasAdapter() throws {
        let packageRoot = try packageRoot()
        let workspaceRoot = packageRoot.deletingLastPathComponent()
        let chemistryRoot = workspaceRoot.appendingPathComponent("AtomLens/Chemistry", isDirectory: true)

        guard FileManager.default.fileExists(atPath: chemistryRoot.path) else {
            throw XCTSkip("AtomLens host app is not present in this checkout.")
        }

        let expectedAdapterRelativePath = "CDKPortAliases.swift"
        var swiftFiles: [String] = []
        let enumerator = FileManager.default.enumerator(at: chemistryRoot, includingPropertiesForKeys: nil)
        while let next = enumerator?.nextObject() as? URL {
            guard next.pathExtension == "swift" else { continue }
            let relative = next.path.replacingOccurrences(of: chemistryRoot.path + "/", with: "")
            swiftFiles.append(relative)
        }

        swiftFiles.sort()
        let unexpected = swiftFiles.filter { $0 != expectedAdapterRelativePath }
        XCTAssertTrue(unexpected.isEmpty, """
        Found AtomLens chemistry source files outside the alias adapter:
        \(unexpected.joined(separator: "\n"))
        Move CDK-derived implementation into CDKSwiftNativePort.
        """)
    }

    private func packageRoot() throws -> URL {
        var url = URL(fileURLWithPath: #filePath)
        // #filePath -> .../Tests/CDKSwiftNativePortTests/PackageBoundaryTests.swift
        url.deleteLastPathComponent()
        url.deleteLastPathComponent()
        url.deleteLastPathComponent()
        guard FileManager.default.fileExists(atPath: url.appendingPathComponent("Package.swift").path) else {
            throw NSError(domain: "PackageBoundaryTests",
                          code: 1,
                          userInfo: [NSLocalizedDescriptionKey: "Could not locate package root from \(#filePath)"])
        }
        return url
    }
}
