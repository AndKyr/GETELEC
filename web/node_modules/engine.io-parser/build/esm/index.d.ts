import encodePacket from "engine.io-parser/build/esm/encodePacket.js";
import decodePacket from "engine.io-parser/build/esm/decodePacket.js";
import { Packet, PacketType, RawData, BinaryType } from "engine.io-parser/build/esm/commons.js";
declare const encodePayload: (packets: Packet[], callback: (encodedPayload: string) => void) => void;
declare const decodePayload: (encodedPayload: string, binaryType?: BinaryType) => Packet[];
export declare const protocol = 4;
export { encodePacket, encodePayload, decodePacket, decodePayload, Packet, PacketType, RawData, BinaryType };
