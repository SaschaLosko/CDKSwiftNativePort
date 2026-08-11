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
            "import SwiftUI",
            "import CoreSpotlight",
            "import QuickLook",
            "import QuickLookThumbnailing"
        ]

        var violations: [String] = []
        let enumerator = FileManager.default.enumerator(at: sourcesRoot,
                                                        includingPropertiesForKeys: nil)

        while let next = enumerator?.nextObject() as? URL {
            guard next.pathExtension == "swift" else { continue }
            let content = try loadUTF8Source(from: next)
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

    func testPackageManifestDeclaresIOSSupport() throws {
        let root = try packageRoot()
        let manifestURL = root.appendingPathComponent("Package.swift")
        let manifest = try String(contentsOf: manifestURL, encoding: .utf8)

        XCTAssertTrue(manifest.contains(".iOS("),
                      "Expected CDKSwiftNativePort Package.swift to declare iOS platform support.")
    }

    func testPackageManifestDeclaresDynamicLibraryProduct() throws {
        let root = try packageRoot()
        let manifestURL = root.appendingPathComponent("Package.swift")
        let manifest = try String(contentsOf: manifestURL, encoding: .utf8)

        XCTAssertTrue(
            manifest.contains("type: .dynamic"),
            "Expected CDKSwiftNativePort to remain an explicitly dynamic library product."
        )
    }

    func testPackageManifestDoesNotCarryAtomLensCrossPlatformCouplingGuidance() throws {
        let root = try packageRoot()
        let manifestURL = root.appendingPathComponent("Package.swift")
        let manifest = try String(contentsOf: manifestURL, encoding: .utf8)

        XCTAssertFalse(manifest.contains("AtomLens app"),
                       "Package manifest should stay package-scoped and not describe host-app coupling.")
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

    func testHostXcodeProjectUsesPackageProductInsteadOfPackageSourceFiles() throws {
        let packageRoot = try packageRoot()
        let workspaceRoot = packageRoot.deletingLastPathComponent()
        let projectURL = workspaceRoot.appendingPathComponent("AtomLens.xcodeproj/project.pbxproj", isDirectory: false)

        guard FileManager.default.fileExists(atPath: projectURL.path) else {
            throw XCTSkip("AtomLens Xcode project is not present in this checkout.")
        }

        let project = try String(contentsOf: projectURL, encoding: .utf8)

        let packageReferenceMarkers = [
            "XCLocalSwiftPackageReference \"CDKSwiftNativePort\"",
            "XCRemoteSwiftPackageReference \"CDKSwiftNativePort\""
        ]
        XCTAssertTrue(packageReferenceMarkers.contains(where: project.contains),
                      "Expected AtomLens to consume the CDKSwiftNativePort Swift package product.")

        let forbiddenMarkers = [
            "CDKSwiftNativePort/Sources/",
            "CDKSwiftNativePort/Tests/"
        ]

        let violations = forbiddenMarkers.filter { project.contains($0) }
        XCTAssertTrue(violations.isEmpty, """
        Found direct Xcode project references to package implementation files:
        \(violations.joined(separator: "\n"))
        Host targets should link the package product, not compile package sources directly.
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

    private func loadUTF8Source(from url: URL) throws -> String {
        let data = try Data(contentsOf: url)
        return String(decoding: data, as: UTF8.self)
    }
}
