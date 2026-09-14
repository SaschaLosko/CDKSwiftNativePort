import Foundation

// CDX tags and encodings follow Revvity's public CDX/CDXML specification:
// https://github.com/Glysade/chemdraw/tree/v1.0.14/docs
// A single object tree backs both encodings, including reaction references.
struct ChemDrawProperty {
    let tag: UInt16
    let name: String
    let text: String
    let bytes: Data

    static func integer<T: FixedWidthInteger>(
        _ tag: UInt16, _ name: String, _ value: T,
        text: String? = nil
    ) -> Self {
        Self(tag: tag, name: name, text: text ?? String(value), bytes: littleEndian(value))
    }

    static func string(_ tag: UInt16, _ name: String, _ value: String) -> Self {
        // One style run explicitly selects the UTF-8 font, including non-ASCII names.
        var bytes = Data()
        for value: UInt16 in [1, 0, 3, 0, 240, 0] { bytes.append(littleEndian(value)) }
        bytes.append(contentsOf: value.utf8)
        return Self(tag: tag, name: name, text: value, bytes: bytes)
    }

    static func references(_ tag: UInt16, _ name: String, _ ids: [UInt32]) -> Self {
        Self(
            tag: tag, name: name, text: ids.map(String.init).joined(separator: " "),
            bytes: ids.reduce(into: Data()) { $0.append(littleEndian($1)) })
    }

    static func coordinates(
        _ tag: UInt16, _ name: String, _ values: [Double],
        binaryOrder: [Int]
    ) throws -> Self {
        var bytes = Data()
        for index in binaryOrder {
            let fixed = (values[index] * 65536).rounded()
            guard fixed.isFinite, fixed >= Double(Int32.min), fixed <= Double(Int32.max) else {
                throw ChemError.unsupported("ChemDraw coordinates exceed the supported page size.")
            }
            bytes.append(littleEndian(Int32(fixed)))
        }
        return Self(tag: tag, name: name, text: values.map { String($0) }.joined(separator: " "), bytes: bytes)
    }
}

struct ChemDrawObject {
    let tag: UInt16
    let name: String
    let id: UInt32
    var properties: [ChemDrawProperty] = []
    var children: [ChemDrawObject] = []

    func xml() -> String {
        let attributes = properties.filter { $0.tag != 0x0700 && $0.tag != 0x0100 }
            .map { " \($0.name)=\"\(Self.escaped($0.text))\"" }.joined()
        var result = "<\(name) id=\"\(id)\"\(attributes)>"
        if tag == 0x8000 {
            result += "<fonttable><font id=\"3\" charset=\"utf-8\" name=\"Arial\"/></fonttable>"
        }
        for property in properties where property.tag == 0x0700 {
            result += "<s font=\"3\" size=\"12\">\(Self.escaped(property.text))</s>"
        }
        result += children.map { $0.xml() }.joined()
        return result + "</\(name)>\n"
    }

    func appendBinary(to data: inout Data) throws {
        data.append(littleEndian(tag == 0x8000 ? UInt16(0) : tag))
        data.append(littleEndian(tag == 0x8000 ? UInt32(0) : id))
        for property in properties {
            guard property.bytes.count <= Int(UInt32.max) else {
                throw ChemError.unsupported("ChemDraw property exceeds the CDX size limit.")
            }
            data.append(littleEndian(property.tag))
            if property.bytes.count >= Int(UInt16.max) {
                data.append(littleEndian(UInt16.max))
                data.append(littleEndian(UInt32(property.bytes.count)))
            } else {
                data.append(littleEndian(UInt16(property.bytes.count)))
            }
            data.append(property.bytes)
        }
        for child in children { try child.appendBinary(to: &data) }
        data.append(littleEndian(UInt16(0)))
    }

    private static func escaped(_ value: String) -> String {
        value.replacingOccurrences(of: "&", with: "&amp;")
            .replacingOccurrences(of: "<", with: "&lt;")
            .replacingOccurrences(of: ">", with: "&gt;")
            .replacingOccurrences(of: "\"", with: "&quot;")
            .replacingOccurrences(of: "\r", with: "&#13;")
            .replacingOccurrences(of: "\n", with: "&#10;")
            .replacingOccurrences(of: "\t", with: "&#9;")
    }
}

func littleEndian<T: FixedWidthInteger>(_ value: T) -> Data {
    var value = value.littleEndian
    return withUnsafeBytes(of: &value) { Data($0) }
}
